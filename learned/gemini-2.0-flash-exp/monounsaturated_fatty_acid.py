"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a MUFA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Find the carbon atoms connected to the carboxyl group
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not matches:
        return False, "No carboxylic acid group found"
    
    carboxyl_atom_idx = matches[0][0] #Get first atom of first match
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_atom_idx)

    
    # Get the carbon connected to carboxyl group (first neighbor not oxygen)
    first_chain_carbon = None
    for neighbor in carboxyl_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            first_chain_carbon = neighbor
            break

    if first_chain_carbon is None:
        return False, "No carbon chain found connected to carboxyl group"

    # Initialize the chain with the first carbon connected to carboxyl group
    chain_atoms = [first_chain_carbon.GetIdx()]
    visited_atoms = {first_chain_carbon.GetIdx()}
    
    
    # Now look for chain of carbons
    stack = [first_chain_carbon.GetIdx()]

    while stack:
        current_idx = stack.pop()
        current_atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in current_atom.GetNeighbors():
             if neighbor.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(current_idx, neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE and neighbor.GetIdx() not in visited_atoms:
                 chain_atoms.append(neighbor.GetIdx())
                 visited_atoms.add(neighbor.GetIdx())
                 stack.append(neighbor.GetIdx())
                 
    chain_mol = Chem.PathToSubmol(mol, chain_atoms)

    # Count double and triple bonds within the chain
    num_double_bonds = len(chain_mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    num_triple_bonds = len(chain_mol.GetSubstructMatches(Chem.MolFromSmarts("C#C")))
    total_unsaturations = num_double_bonds + num_triple_bonds

    if total_unsaturations != 1:
        return False, f"Molecule has {total_unsaturations} double/triple bonds in the fatty acid chain, should have exactly 1"

    # Count the number of carbons in the chain
    c_count = len(chain_atoms)

    # Check for fatty acid chain length (at least 4 carbons in the chain, as smallest MUFA is butenoic acid)
    if c_count < 4:
        return False, f"Carbon chain length ({c_count}) is too short for a fatty acid"


    return True, "Monounsaturated fatty acid identified"