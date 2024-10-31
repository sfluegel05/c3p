from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_spirostanyl_glycoside(smiles: str):
    """
    Determines if a molecule is a spirostanyl glycoside based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a spirostanyl glycoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS pattern for spirostan core with more flexible ring pattern
    spirostan_pattern = Chem.MolFromSmarts('[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]2~[#6]~[#6]4~[#6]3~[#6]~[#6]5~[#6]4O[#6]6[#6][#6][#6][#6]O6~[#6]5~[#6]1')

    # More flexible SMARTS pattern for glycoside moiety
    glycoside_pattern = Chem.MolFromSmarts('O[CH]1O[CH2,CH][CH]([OH,O])[CH]([OH,O])[CH]1[OH,O]')

    # Check for spirostan core
    if not mol.HasSubstructMatch(spirostan_pattern):
        return False, "No spirostan core found"

    # Check for glycoside
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycoside moiety found"

    # Check for O-glycosidic linkage between spirostan and sugar
    glycosidic_link = Chem.MolFromSmarts('[#6]~[#6]~O~[#6]1O[#6][#6]([#6]1)[OH,O]')
    if not mol.HasSubstructMatch(glycosidic_link):
        return False, "No O-glycosidic linkage found"

    # Count rings
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    
    # Need at least 6 rings (4 steroid + 2 spiro)
    if len(rings) < 6:
        return False, "Insufficient number of rings"

    # Count oxygen atoms for glycoside verification
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    if o_count < 5:
        return False, "Insufficient oxygen atoms for glycoside"

    return True, "Valid spirostanyl glycoside structure found"
# Pr=None
# Recall=0.0