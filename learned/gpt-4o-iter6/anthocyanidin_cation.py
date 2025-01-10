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
    
    # Basic flavylium core pattern, with more flexibility in substituents
    flavylium_pattern = Chem.MolFromSmarts("[o+]1c2ccccc2c(O)c2cc(O)ccc12")  
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "Does not match the anthocyanidin cation core (flexible flavylium) pattern"
    
    # Ensure positive charge is found - common in these cations
    if not any(atom.GetFormalCharge() == 1 for atom in mol.GetAtoms()):
        return False, "Missing positive charge expected in anthocyanidin cations"
    
    # Check for a reasonable count of oxygen atoms, as these are often oxygenated
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Insufficient number of oxygen atoms ({o_count}) for anthocyanidin cation"

    # Check for some common glycosidic linkages, though not exclusively defining
    glycoside_links = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O")
    if mol.HasSubstructMatch(glycoside_links):
        return True, "Contains anthocyanidin cation with glycosidic linkage"

    return True, "Contains anthocyanidin cation core with sufficient oxygenation"