"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    Oligosaccharides generally contain 3-10 monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for pyranose/furanose rings (considering common structural variety)
    pyranose_pattern = Chem.MolFromSmarts("[C&R]1OC(C(C1)O)O")
    furanose_pattern = Chem.MolFromSmarts("[C&R]1O[C&R]([C&R]1)O")
    
    pyranose_matches = len(mol.GetSubstructMatches(pyranose_pattern))
    furanose_matches = len(mol.GetSubstructMatches(furanose_pattern))
    total_ring_matches = pyranose_matches + furanose_matches
    
    # Improved glycosidic bond pattern
    glycosidic_pattern = Chem.MolFromSmarts("[C&R]1O[C&R][OX2&R]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    # Count the number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Adjust criteria based on expected characteristics of an oligosaccharide
    if 2 <= total_ring_matches <= 10 and 10 <= o_count <= 50 and glycosidic_matches >= 2:
        return True, "Contains multiple saccharide-like structures with glycosidic linkages"
    else:
        reason = (f"Characterization doesn't match oligosaccharide: "
                  f"{total_ring_matches} rings, {o_count} oxygens, {glycosidic_matches} glycosidic-like bonds")
        return False, reason