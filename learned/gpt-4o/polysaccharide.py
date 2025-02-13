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

    # Define a SMARTS pattern that matches a common pyranose ring
    pyranose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)CCO1")
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)

    if len(pyranose_matches) <= 10:
        return False, f"Found {len(pyranose_matches)} pyranose rings, need more than 10"

    # Define a SMARTS pattern for glycosidic linkage
    glycosidic_pattern = Chem.MolFromSmarts("C1OC(CCO1)C(O)C")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    # Ensure that both pyranose rings and linkages exceed the thresholds
    if len(glycosidic_matches) < 10:
        return False, f"Found {len(glycosidic_matches)} potential glycosidic linkages, need more than 10"

    return True, "Contains more than ten pyranose rings linked via glycosidic bonds"