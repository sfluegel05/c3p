"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:24169 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound in which two monosaccharides are joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify monosaccharide rings
    monosaccharide_rings = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        oxygen_count = sum(atom.GetAtomicNum() == 8 for atom in ring_atoms)
        hydroxyl_count = sum(atom.GetTotalNumHs(onlyExplicit=False) == 1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        anomeric_carbons = [atom.GetIdx() for atom in ring_atoms if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP2]
        
        if len(ring) in [5, 6] and oxygen_count >= 1 and hydroxyl_count >= 1 and len(anomeric_carbons) == 1:
            monosaccharide_rings.append(ring)
    
    # Check for glycosidic bond
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CH1X4][CH1X4][OX2]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    
    if len(monosaccharide_rings) != 2 or len(glycosidic_bond_matches) == 0:
        return False, "Does not contain two monosaccharide rings connected by a glycosidic bond"
    
    # Verify the glycosidic bond connection
    for match in glycosidic_bond_matches:
        ring1, ring2 = None, None
        for ring in monosaccharide_rings:
            if match[1] in ring:
                ring1 = ring
            if match[2] in ring:
                ring2 = ring
        
        if ring1 and ring2:
            break
    else:
        return False, "Glycosidic bond does not connect monosaccharide rings"
    
    return True, "Contains two monosaccharide rings connected by a glycosidic bond"