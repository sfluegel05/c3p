from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthracycline_antibiotic(smiles: str):
    """
    Determines if a molecule is an anthracycline antibiotic based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anthracycline antibiotic, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for tetrahydronaphthacenedione core structure
    core_pattern = Chem.MolFromSmarts('O=C1C=C2C(=O)c3c(O)cccc3C(=O)c2c(O)c1')

    # SMARTS pattern for amino sugar (daunosamine)
    sugar_pattern = Chem.MolFromSmarts('[CH3][CH]1O[CH2][CH]([NH2])[CH]([OH])[CH]1[OH]')

    # SMARTS pattern for glycosidic linkage
    glycosidic_pattern = Chem.MolFromSmarts('[#6]-[#8]-[CH]1-[CH2]-[CH]([NH2])-[CH]-[CH]-[#8]-1')

    # Check for required structural features
    has_core = mol.HasSubstructMatch(core_pattern)
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    has_glycosidic = mol.HasSubstructMatch(glycosidic_pattern)

    # Ring count to verify tetracyclic system
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()

    if not has_core and ring_count < 4:
        return False, "Missing tetracyclic ring system with quinone structure"
    
    if not (has_sugar or has_glycosidic):
        return False, "Missing daunosamine sugar moiety or glycosidic linkage"

    # Count oxygen atoms to ensure presence of hydroxyl groups
    oxygen_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if oxygen_count < 5:
        return False, "Insufficient oxygen-containing groups"

    # Additional check for aromatic rings
    aromatic_rings = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_rings < 2:
        return False, "Insufficient aromatic system"

    return True, "Contains characteristic anthracycline antibiotic structural features"
# Pr=1.0
# Recall=0.6111111111111112