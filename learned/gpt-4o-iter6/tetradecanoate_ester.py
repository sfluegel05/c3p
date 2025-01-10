"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

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

    # Define a SMARTS pattern for a tetradecanoate ester
    # Looking for general ester pattern with 14-carbon chain
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)OC")  # 14-carbon chain including ester formation
    
    if mol.HasSubstructMatch(tetradecanoate_pattern):
        # Get matches
        matches = mol.GetSubstructMatches(tetradecanoate_pattern)
        for match in matches:
            c_count = sum(1 for atom_idx in match if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
            if c_count == 14:  # Validating the exact carbon count expected
                return True, "Valid tetradecanoate ester group found"
        
    return False, "No valid tetradecanoate ester group found or incorrect carbon count"