"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length of less than C6
    (maximum of 5 carbons in the alkyl chain) and no non-hydrocarbon substituents other than
    the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group and extract the chain
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CH2;!$(C=O)]-[C](=O)[OH]")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not matches:
        carboxylic_acid_pattern = Chem.MolFromSmarts("[CH;!$(C=O)]-[C](=O)[OH]") # Check also for -CH- if it is not a terminal carbon
        matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
        if not matches:
          return False, "No carboxylic acid group found"

    alpha_carbon_idx = matches[0][0]
    
    def get_chain_atoms(mol, start_atom_idx, visited_atoms):
        """
        Recursively finds the carbon atoms in the alkyl chain connected to the carboxylic acid group.
        Only single bonds are followed.
        """
        visited_atoms.add(start_atom_idx)
        chain_atoms = [start_atom_idx]

        neighbors = []
        for neighbor in mol.GetAtomWithIdx(start_atom_idx).GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(start_atom_idx, neighbor_idx)

            if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited_atoms and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
              neighbors.append(neighbor_idx)

        for neighbor_idx in neighbors:
          chain_atoms.extend(get_chain_atoms(mol, neighbor_idx, set(visited_atoms)))


        return chain_atoms

    chain_atoms = get_chain_atoms(mol, alpha_carbon_idx,set())
    
    #check if there are any non-hydrocarbons in the chain atoms. Oxygen is OK if part of the -OH group
    for atom_idx in chain_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1 and not (atom.GetAtomicNum() == 8 and atom.GetIdx() != mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O"))[0][2] ):
              return False, "Non-hydrocarbon substituent found on the chain"

    # Calculate chain length
    chain_length = len(chain_atoms)
    
    if chain_length > 5 :
        return False, "More than 5 carbons in the alkyl chain"
    
    
    return True, "Short-chain fatty acid criteria met"