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
    
    # Look for pyranose (6-membered) and furanose (5-membered) sugar rings
    pyranose_pattern = Chem.MolFromSmarts("C1OC[C@H](O)[C@@H](O)C1")
    furanose_pattern = Chem.MolFromSmarts("C1OC[C@H](O)C1")
    
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    total_ring_matches = len(pyranose_matches) + len(furanose_matches)
    
    # Improve glycosidic bond pattern (O-C-O linkage)
    glycosidic_pattern = Chem.MolFromSmarts("C-O-C")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    # Count the number of simple ether and hydroxyl groups - relevant to sugar structures
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Determine if the molecule is likely an oligosaccharide
    if 3 <= total_ring_matches <= 10 and 6 <= o_count and len(glycosidic_matches) >= 2:
        return True, "Contains multiple saccharide-like structures with glycosidic linkages"
    else:
        reason = (f"Characterization doesn't match oligosaccharide: "
                  f"{total_ring_matches} rings, {o_count} oxygens, {len(glycosidic_matches)} glycosidic-like bonds")
        return False, reason