"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain C18 fatty acid with two C=C double bonds and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for exactly two C=C double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 2:
         return False, f"Found {len(double_bond_matches)} double bonds, need exactly 2"

    # Check for a long straight chain of 18 C atoms, where the double bonds must reside
    # Find the carbon chain atoms first, then verify the chain length
    all_carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]

    # Get the carbon atoms connected to the carboxylic acid
    carboxylic_carbons = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O[C]"))
    if not carboxylic_carbons:
       return False, "No carbon chain directly attached to the carboxylic acid"
    carboxylic_carbon_atom = mol.GetAtomWithIdx(carboxylic_carbons[0][0])

    # Find the longest chain by following the carbons from the acid end
    carbons_in_chain = [carboxylic_carbon_atom]
    current_atom = carboxylic_carbon_atom
    previous_atom = None
    chain_length = 1

    while chain_length < 18:
        # Get the neighbors of the current atom, except the previous atom
        neighbors = [neighbor for neighbor in current_atom.GetNeighbors()
                     if neighbor.GetAtomicNum() == 6 and neighbor != previous_atom]
        
        # If there is exactly one neighbor in the chain, then follow it
        if len(neighbors) == 1:
            next_atom = neighbors[0]
            carbons_in_chain.append(next_atom)
            previous_atom = current_atom
            current_atom = next_atom
            chain_length += 1
        else:
            # If it doesn't have exactly one carbon neighbor, this is not a straight chain
            return False, "Not a straight chain with 18 carbons"

    #verify that all the double bonds are on the 18-carbon chain
    double_bonds_on_chain = 0
    for match in double_bond_matches:
        for atom_idx in match:
            if mol.GetAtomWithIdx(atom_idx) in carbons_in_chain:
                double_bonds_on_chain +=1
                break
    if double_bonds_on_chain != 4:
      return False, "Double bonds not on the 18 carbon chain"

    # Check for non-cyclic structure
    if mol.GetRingInfo().NumRings() > 0:
         return False, "Molecule contains a ring structure"
    

    return True, "Contains 18 carbons in a straight chain, a carboxylic acid group and two double bonds"