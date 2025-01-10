"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem.rdchem import Atom

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

    # Count the number of saccharide rings using pyranose (C1CCO1) or furanose (C1CO1) patterns
    pyranose_pattern = Chem.MolFromSmarts("[C&R]1O[C&R][C&R][C&R][C&R]O1")
    furanose_pattern = Chem.MolFromSmarts("[C&R]1O[C&R][C&R]O1")
    
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    total_ring_matches = len(pyranose_matches) + len(furanose_matches)
    
    # Count the number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Count the number of glycosidic linkages (C-O-C patterns)
    glycosidic_pattern = Chem.MolFromSmarts("[C&R]O[C&R]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    # Criteria for oligosaccharide classification
    if total_ring_matches >= 2 and 10 <= o_count <= 40 and len(glycosidic_matches) >= 1:
        return True, "Contains features characteristic of an oligosaccharide (multiple saccharide rings and glycosidic bonds)"
    else:
        reason = (f"Characterization doesn't match oligosaccharide: "
                  f"{total_ring_matches} rings, {o_count} oxygens, {len(glycosidic_matches)} glycosidic-like bonds")
        return False, reason