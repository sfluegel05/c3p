"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is defined as a compound with multiple monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pyranose and furanose SMARTS pattern
    pyranose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(OC(O)C1)CO")
    furanose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)OC1")

    # Check for multiple monosaccharide units
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    total_units = len(pyranose_matches) + len(furanose_matches)

    if total_units < 2:
        return False, f"Insufficient monosaccharide units, found {total_units}"

    # Check for glycosidic linkages (ether linkages between monosaccharides)
    ether_pattern = Chem.MolFromSmarts("OCCO")
    ether_matches = mol.GetSubstructMatches(ether_pattern)

    if len(ether_matches) < total_units - 1:
        return False, f"Insufficient glycosidic linkages for an oligosaccharide, found {len(ether_matches)}"

    return True, "Molecule has multiple monosaccharide units connected by glycosidic linkages"