"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are oxygenated derivatives of flavylium ions.

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
    
    # Adjusted SMARTS pattern for anthocyanidin core with more flexibility in oxygen substitutions
    anthocyanidin_core_pattern = Chem.MolFromSmarts("[o+]1c2cc(O)c(O)c(O)c2cc3c(O)ccc(O)c13")
    if not mol.HasSubstructMatch(anthocyanidin_core_pattern):
        return False, "Does not match the flexible anthocyanidin cation core pattern"
    
    # Verify the presence of a positive charge on an oxygen atom
    charge_found = any(atom.GetFormalCharge() == 1 and atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not charge_found:
        return False, "No positive charge found in expected anthocyanidin cation"
    
    # Count the number of oxygens and confirm it includes at least three
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Insufficient number of oxygen atoms ({o_count}) for anthocyanidin cation"

    # Identification of potential glycosidic linkages – flexible pattern to accommodate different sugars
    glycoside_pattern = Chem.MolFromSmarts("O[C@H1]C(O)C(O)C[*]")
    if mol.HasSubstructMatch(glycoside_pattern):
        return True, "Contains anthocyanidin cation with glycosidic linkage"
    
    return True, "Contains anthocyanidin cation core with oxygenated structure and charge characteristics"