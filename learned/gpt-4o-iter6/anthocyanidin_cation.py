"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are oxygenated derivatives of flavylium (2-phenylchromenylium).

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
    
    # Look for 2-phenylchromenylium core with positive charge on oxygen or carbon
    flavylium_pattern = Chem.MolFromSmarts("[O+]1[cH]cc2c(O)ccc(O)c2c1")  
    # The pattern captures the positively charged oxygen(O+) or carbon in the flavylium cation core
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "Does not match the anthocyanidin cation core (flavylium) pattern"
    
    # Check for sufficient oxygenation, typical for anthocyanidins
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient number of oxygen substituents for anthocyanidin cation"

    return True, "Contains 2-phenylchromenylium core with sufficient oxygenation"