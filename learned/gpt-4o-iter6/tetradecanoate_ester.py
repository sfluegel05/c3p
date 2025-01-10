"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is an ester formed from tetradecanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a feltative comprehensive SMARTS pattern to match tetradecanoate ester.
    # The goal is to match C(=O)O after a 13-carbon chain: "CCCCCCCCCCCCC(=O)O"
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)O") 

    # Check for the presence of the tetradecanoate ester structure
    if not mol.HasSubstructMatch(tetradecanoate_pattern):
        return False, "No tetradecanoate ester group found"

    # Validate by the number of carbon atoms in the tetradecanoate portion
    for match in mol.GetSubstructMatches(tetradecanoate_pattern):
        chain_count = sum(1 for atom_idx in match if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
        if chain_count != 14:  # 13 from chain, 1 from the carbonyl C
            return False, f"Incorrect carbon count ({chain_count}) for tetradecanoate ester structure, expected 14"

    return True, "Valid tetradecanoate ester functionality identified"