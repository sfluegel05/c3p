"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: long-chain fatty acid (C13 to C22)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13 to C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)[OX1H0,OX2H1]") # check for -COOH and -COO-
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Get the carbon atoms in the carboxylic acid
    carboxylic_carbon = mol.GetSubstructMatch(acid_pattern)[0]
    carboxylic_atom = mol.GetAtomWithIdx(carboxylic_carbon)
    
    # Get adjacent carbon atoms
    chain_carbons = []
    for neighbor in carboxylic_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
          chain_carbons.append(neighbor.GetIdx())

    # Recursive function to collect all connected carbon atoms
    def find_connected_carbons(mol, start_idx, visited=None):
        if visited is None:
            visited = set()
        if start_idx in visited:
            return visited

        visited.add(start_idx)
        
        current_atom = mol.GetAtomWithIdx(start_idx)
        for neighbor in current_atom.GetNeighbors():
          if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carboxylic_carbon: # Ignore the carboxylic acid carbon
            find_connected_carbons(mol, neighbor.GetIdx(), visited)
        return visited

    all_chain_carbons = set()
    for carbon_idx in chain_carbons:
      all_chain_carbons.update(find_connected_carbons(mol, carbon_idx))


    # Count the number of carbons
    carbon_count = len(all_chain_carbons)


    # Check if the chain length is within the range of 13 to 22 carbons.
    if carbon_count < 13 or carbon_count > 22:
        return False, f"Chain length is {carbon_count}, not between 13 and 22 carbons"

    return True, f"Molecule is a long-chain fatty acid with a chain length of {carbon_count} carbons."