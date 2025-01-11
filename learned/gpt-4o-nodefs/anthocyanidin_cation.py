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

    # Look for flavylium cation core structure
    flavylium_pattern = Chem.MolFromSmarts("c1(c2ccc(O)cc2)C=[O+]c2cc(O)cc(O)c12")
    if mol.HasSubstructMatch(flavylium_pattern):
        return True, "Contains flavylium cation structure characteristic of anthocyanidin cations"
    else:
        return False, "Does not contain flavylium cation structure"