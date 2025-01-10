"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: CHEBI:24474 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic chromene (benzopyran) core with aryl at position 3
    # The =O allows for ketone at position 4 (common in isoflavones)
    isoflavonoid_pattern = Chem.MolFromSmarts('[cR1]1[cR1][cR1][cR1][cR1][cR1]1-[CR1]2=C([ar5])[CR1](=O)[cR1]1[cR1][cR1][cR1][cR1][cR1]1O2')
    
    # Alternative pattern for isoflavanones (saturated C=O bond)
    isoflavanone_pattern = Chem.MolFromSmarts('[cR1]1[cR1][cR1][cR1][cR1][cR1]1-[CR1]2[CR1]([ar5])[CR1](=O)[cR1]1[cR1][cR1][cR1][cR1][cR1]1O2')
    
    if isoflavonoid_pattern is None or isoflavanone_pattern is None:
        return None, "Error in SMARTS patterns"

    # Check for required substructures
    has_isoflavone = mol.HasSubstructMatch(isoflavonoid_pattern)
    has_isoflavanone = mol.HasSubstructMatch(isoflavanone_pattern)
    
    if not (has_isoflavone or has_isoflavanone):
        return False, "No isoflavonoid core structure found"

    # Count rings to ensure we have at least 3 (two from chromene, one from aryl)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring count for isoflavonoid structure"

    # Most isoflavonoids have multiple oxygen-containing substituents
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Insufficient oxygen atoms for typical isoflavonoid"

    # Count aromatic rings to ensure we have the benzopyran and at least one aryl
    aromatic_rings = len(mol.GetSubstructMatches(Chem.MolFromSmarts('a1aaaaa1')))
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings for isoflavonoid"

    if has_isoflavone:
        return True, "Contains isoflavone core structure with aryl substituent at position 3"
    else:
        return True, "Contains isoflavanone core structure with aryl substituent at position 3"