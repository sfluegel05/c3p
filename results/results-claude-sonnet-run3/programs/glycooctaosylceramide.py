from rdkit import Chem
from rdkit.Chem import AllChem
import re

def is_glycooctaosylceramide(smiles: str):
    """
    Determines if a molecule is a glycooctaosylceramide based on SMILES string.
    A glycooctaosylceramide consists of a glycooctaosyl moiety attached to ceramide.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glycooctaosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Check for valid SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ceramide core structure
    # Looking for long alkyl chain, sphingosine base, and amide linkage
    ceramide_pattern = r"CCCCCCCCCCCCC\\C=C\\[C@@H]\(O\)[C@H].*NC\(=O\)"
    if not re.search(ceramide_pattern, smiles):
        return False, "Missing ceramide core structure"

    # Count sugar rings
    rings = mol.GetRingInfo()
    sugar_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    
    if len(sugar_rings) != 8:
        return False, f"Found {len(sugar_rings)} sugar rings, need exactly 8"

    # Check for glycosidic linkages between sugars
    glycosidic_pattern = r"O[C@H]|O[C@@H]"
    glycosidic_count = len(re.findall(glycosidic_pattern, smiles))
    if glycosidic_count < 14:  # Need sufficient linkages for 8 connected sugars
        return False, "Insufficient glycosidic linkages"

    # Verify presence of hydroxyl groups on sugars
    hydroxyl_pattern = r"[C@H]\(O\)|[C@@H]\(O\)"
    hydroxyl_count = len(re.findall(hydroxyl_pattern, smiles))
    if hydroxyl_count < 16:  # Each sugar should have multiple hydroxyls
        return False, "Insufficient hydroxyl groups"

    # Check for characteristic sugar CH2OH groups
    ch2oh_pattern = r"CO\)"
    ch2oh_count = len(re.findall(ch2oh_pattern, smiles))
    if ch2oh_count < 6:  # Most sugars should have CH2OH
        return False, "Insufficient CH2OH groups"

    # Check for N-acetyl groups
    nacetyl_pattern = r"NC\(=O\)C"
    nacetyl_count = len(re.findall(nacetyl_pattern, smiles))
    if nacetyl_count < 1:
        return False, "Missing N-acetyl groups"

    # Check for the ceramide-sugar linkage
    if not re.search(r"CO[C@@H].*[C@H].*NC\(=O\)", smiles):
        return False, "Missing proper ceramide-sugar linkage"

    # If all checks pass
    return True, "Structure contains ceramide core with 8 connected sugar units and appropriate linkages"
# Pr=None
# Recall=0.0