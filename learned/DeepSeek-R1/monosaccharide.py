"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:35724 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde/ketone with ≥3 carbons or a cyclic hemiacetal, without glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for at least three carbons
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 3:
        return False, "Less than three carbons"
    
    # Count hydroxyl groups (O with H)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0)
    if hydroxyl_count < 2:
        return False, f"Only {hydroxyl_count} hydroxyl groups"
    
    # Check for aldehyde groups [CX3H1]=O
    has_aldehyde = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3H1]=O'))) > 0
    
    # Check for ketones not adjacent to O (exclude acids/esters)
    has_ketone = False
    ketone_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[#6]'))
    for match in ketone_matches:
        carbon = mol.GetAtomWithIdx(match[0])
        if not any(n.GetAtomicNum() == 8 for n in carbon.GetNeighbors() if n.GetIdx() != match[1]):
            has_ketone = True
            break
    
    # Check for rings with oxygen (potential cyclic hemiacetal)
    rings = mol.GetRingInfo().AtomRings()
    has_ring_oxygen = any(any(mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in ring) for ring in rings)
    
    if not (has_aldehyde or has_ketone or has_ring_oxygen):
        return False, "No aldehyde/ketone or cyclic hemiacetal"
    
    # Check for glycosidic bonds (O connected to two non-ring carbons)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and len(atom.GetNeighbors()) == 2:
            if not atom.IsInRing() and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors()):
                return False, "Glycosidic bond present"
    
    return True, "Meets monosaccharide criteria"