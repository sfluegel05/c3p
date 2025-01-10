"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester contains a tetradecanoate moiety, which is a 14-carbon saturated fatty acid ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern intended to catch -CCCCCCCCCCCCCC(=O)O- (14 carbons ester link)
    # This allows for flexibility in downstream esterification (multi-ester structures etc.)
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")

    # Check for presence of the tetradecanoate moiety
    if mol.HasSubstructMatch(tetradecanoate_pattern):
        return True, "Contains tetradecanoate moiety"

    return False, "Does not contain tetradecanoate moiety"