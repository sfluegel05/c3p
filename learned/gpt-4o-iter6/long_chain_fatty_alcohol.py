"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a carbon chain length ranging from C13 to C22 with
    at least one hydroxyl (-OH) group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for hydroxyl group (-OH) presence
    oh_group = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(oh_group):
        return False, "No hydroxyl group found"
    
    # Function to find the longest chain of carbon atoms
    def longest_carbon_chain(mol):
        max_length = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Only consider carbon atoms
                lengths = []
                rdmolops.GetBondsInOrder(atom, molecule=mol, atomfun=lambda x: x.GetAtomicNum() == 6, lengthfun=len, results=lengths)
                if lengths:
                    max_length = max(lengths + [max_length])
        return max_length
    
    # Get the longest carbon chain in the molecule
    max_chain_length = longest_carbon_chain(mol)
    
    # Check if any chain meets the C13 to C22 requirement
    if 13 <= max_chain_length <= 22:
        return True, f"Contains a carbon chain of length {max_chain_length} and a hydroxyl group"
    else:
        return False, f"Carbon chain length of {max_chain_length}, required between 13 and 22"