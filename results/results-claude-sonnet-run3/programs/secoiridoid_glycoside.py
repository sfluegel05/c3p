from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_secoiridoid_glycoside(smiles: str):
    """
    Determines if a molecule is a secoiridoid glycoside.
    A secoiridoid glycoside is a glycoside of any secoiridoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secoiridoid glycoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of pyranose sugar (glycoside)
    sugar_pattern = Chem.MolFromSmarts("[CH1,CH2]1O[CH1]([CH1]([CH1]([CH1]1)[OH1,OH0])[OH1,OH0])[OH1,OH0]")
    
    # Core secoiridoid patterns:
    # Pattern 1: Cyclopentane with unsaturated lactone or ester
    secoiridoid_core1 = Chem.MolFromSmarts("C1CC(=O)OCC1")
    
    # Pattern 2: Core with opened iridoid structure containing methylester
    secoiridoid_core2 = Chem.MolFromSmarts("C=CC1OC=C(C(=O)OC)C1")
    
    # Pattern 3: Alternative opened iridoid core with ester
    secoiridoid_core3 = Chem.MolFromSmarts("C1OC=CC(CC(=O)O)=C1")
    
    # Pattern 4: Another secoiridoid variation
    secoiridoid_core4 = Chem.MolFromSmarts("C1OC=C(C(=O)O)C1")

    # Check for presence of both sugar and secoiridoid core
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    has_core = (mol.HasSubstructMatch(secoiridoid_core1) or 
                mol.HasSubstructMatch(secoiridoid_core2) or
                mol.HasSubstructMatch(secoiridoid_core3) or
                mol.HasSubstructMatch(secoiridoid_core4))

    # Check glycosidic linkage patterns
    glycosidic_bond1 = Chem.MolFromSmarts("C1OCCC(OC2OCCCC2)C1")
    glycosidic_bond2 = Chem.MolFromSmarts("C1OC=CC(OC2OCCCC2)C1")
    has_linkage = (mol.HasSubstructMatch(glycosidic_bond1) or 
                   mol.HasSubstructMatch(glycosidic_bond2))

    if not has_sugar:
        return False, "No glycoside moiety found"
    
    if not has_core:
        return False, "No secoiridoid core structure found"
        
    if not has_linkage:
        return False, "No glycosidic linkage found"

    # Check molecular weight range
    mol_weight = Descriptors.ExactMolWt(mol)
    if mol_weight < 300 or mol_weight > 1000:
        return False, f"Molecular weight {mol_weight:.1f} outside typical range for secoiridoid glycosides"

    return True, "Contains secoiridoid core with glycoside moiety linked via glycosidic bond"
# Pr=None
# Recall=0.0