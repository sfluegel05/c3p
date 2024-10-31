from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycononaosylceramide(smiles: str):
    """
    Determines if a molecule is a glycononaosylceramide (oligoglycosylceramide with 9 sugar units).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glycononaosylceramide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for ceramide moiety
    ceramide_patterns = [
        # Long chain base with amide
        Chem.MolFromSmarts('CCCCCCCCCCCCC[CH]=C[CH]([CH](CO)NC=O)O'),
        # Alternative ceramide pattern
        Chem.MolFromSmarts('CCCCCCCCCCCCC[CH2][CH]([CH](CO)NC=O)O'),
        # More general pattern
        Chem.MolFromSmarts('CCCCCCCCCC[CH2][CH2][CH]=C[CH]([CH](CO)NC=O)O')
    ]
    
    has_ceramide = False
    for pattern in ceramide_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_ceramide = True
            break
            
    if not has_ceramide:
        return False, "Missing ceramide moiety"

    # Count sugar rings
    sugar_patterns = [
        # Pyranose ring with hydroxyl groups
        Chem.MolFromSmarts('[C]1[C][C]([O,N])[C]([O,N])[C]([O,N])[O]1'),
        # Alternative pyranose pattern
        Chem.MolFromSmarts('[C]1[O][C]([CH2][O,N])[C]([O,N])[C]([O,N])[C]1[O,N]')
    ]

    sugar_rings = set()
    for pattern in sugar_patterns:
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                sugar_rings.update(match)
    
    # Count unique sugar rings by looking at ring atoms
    unique_rings = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:  # Only count 6-membered rings
            ring_atoms = set(ring)
            if any(atom in sugar_rings for atom in ring_atoms):
                unique_rings += 1

    if unique_rings != 9:
        return False, f"Contains {unique_rings} sugar rings instead of 9"

    # Check for glycosidic linkages between sugars
    glycosidic_pattern = Chem.MolFromSmarts('[C;R]-[O;!R]-[C;R]')
    if glycosidic_pattern is None:
        return None, "Invalid glycosidic pattern SMARTS"
        
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    if len(glycosidic_matches) < 8:  # Need at least 8 glycosidic bonds to connect 9 sugars
        return False, f"Contains {len(glycosidic_matches)} glycosidic bonds instead of at least 8"

    return True, "Contains ceramide moiety and 9 sugar units connected by glycosidic bonds"
# Pr=1.0
# Recall=0.3333333333333333