"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    Polysaccharides consist of more than ten monosaccharide units linked by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern for recognizing hexose/pyranose rings
    # This includes stereochemistry and common variations
    pyranose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)O1")
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)

    if len(pyranose_matches) <= 10:
        return False, f"Found {len(pyranose_matches)} pyranose rings, need more than 10"

    # Enhanced pattern for glycosidic linkages
    # Acknowledges variability in C-O-C connectivity with potential stereochemistry
    glycosidic_pattern = Chem.MolFromSmarts("C[O]-[C@]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    if len(glycosidic_matches) < 10:
        return False, f"Found {len(glycosidic_matches)} potential glycosidic linkages, need more than 10"

    return True, "Contains more than ten pyranose rings linked via glycosidic bonds"