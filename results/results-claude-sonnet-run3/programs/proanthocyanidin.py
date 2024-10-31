from rdkit import Chem
from rdkit.Chem import AllChem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin (condensed tannin).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for minimum number of atoms (proanthocyanidins are oligomers)
    if mol.GetNumAtoms() < 30:
        return False, "Too small to be a proanthocyanidin"

    # Count number of phenol groups (OH attached to aromatic ring)
    phenol_pattern = Chem.MolFromSmarts('c[OH]')
    phenol_matches = len(mol.GetSubstructMatches(phenol_pattern))
    if phenol_matches < 4:
        return False, "Insufficient phenol groups"

    # Look for characteristic flavan-3-ol core structure
    flavan_pattern = Chem.MolFromSmarts('O[CH]1Cc2c(O)cc(O)cc2O[CH]1c1ccc(O)cc1')
    if not mol.HasSubstructMatch(flavan_pattern):
        return False, "Missing flavan-3-ol core structure"

    # Check for characteristic C-C and C-O-C linkages between units
    # 4->8 linkage pattern
    c4_c8_pattern = Chem.MolFromSmarts('[CH]1[CH](O)[CH](Oc2cc(O)cc(O)c12)')
    # 4->6 linkage pattern
    c4_c6_pattern = Chem.MolFromSmarts('[CH]1[CH](O)[CH](Oc2c(O)cc(O)cc12)')
    
    has_characteristic_linkage = (mol.HasSubstructMatch(c4_c8_pattern) or 
                                mol.HasSubstructMatch(c4_c6_pattern))
    
    if not has_characteristic_linkage:
        return False, "Missing characteristic 4->8 or 4->6 linkage"

    # Count the number of flavan units (approximate)
    ring_info = mol.GetRingInfo()
    benzene_rings = sum(1 for ring in ring_info.AtomRings() if len(ring) == 6)
    
    # Each flavan unit typically has 2 benzene rings
    estimated_units = benzene_rings // 2
    
    if estimated_units < 2:
        return False, "Insufficient number of flavan units"

    return True, f"Proanthocyanidin with approximately {estimated_units} flavan units"
# Pr=1.0
# Recall=0.8421052631578947