"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is characterized by a 2-phenylchromenylium core structure
    with a positively charged oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised SMARTS pattern for flavylium core of anthocyanidins
    # Includes a positively charged oxygen part of a 2-phenylchromenylium system
    flavylium_pattern = Chem.MolFromSmarts(
        "[O+]1c2cc(O)cc(O)c2[o+]c2c(cccc2)c2ccccc12"
    )
    
    # Check for presence of flavylium cation core structure
    if mol.HasSubstructMatch(flavylium_pattern):
        return True, "Contains flavylium cation structure characteristic of anthocyanidin cations"
    else:
        return False, "Does not contain flavylium cation structure"