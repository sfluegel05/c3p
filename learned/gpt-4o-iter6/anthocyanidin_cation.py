"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidins are oxygenated derivatives of flavylium cations.

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
    
    # Flexible anthocyanidin core pattern with flavylium character
    # Allow aromatic variation and oxygen substitutions
    flavylium_resonance_pattern = Chem.MolFromSmarts("[o+]1c2ccccc2-c(O)cc(O)c3cc(O)ccc13")
    if not mol.HasSubstructMatch(flavylium_resonance_pattern):
        return False, "Does not match the flexible anthocyanidin cation core pattern"
    
    # Verify the presence of a charge, emphasizing positive charge typical for flavylium-derived cations
    if not any(atom.GetFormalCharge() == 1 for atom in mol.GetAtoms()):
        return False, "No positive charge found in expected anthocyanidin cation"
    
    # Count the number of oxygens to confirm oxygenation
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Insufficient number of oxygen atoms ({o_count}) for anthocyanidin cation"

    # General catch for any acceptable glycosidic substructures without strict stereochemistry
    # Flexible glycoside patterns 
    glycoside_pattern = Chem.MolFromSmarts("O[*]C1C(O)C(O)C[C@H]1") 
    if mol.HasSubstructMatch(glycoside_pattern):
        return True, "Contains anthocyanidin cation with flexible glycosidic linkage"
    
    return True, "Contains anthocyanidin cation core with proper oxygenation and charge characteristics"