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

    # Define more inclusive pyranose and furanose SMARTS patterns for sugars
    pyranose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)")
    furanose_pattern = Chem.MolFromSmarts("O1[C@@H]([C@@H](O)[C@H](O)[C@H]1O)")

    # Match monosaccharide units (pyranose and/or furanose)
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    total_units = len(pyranose_matches) + len(furanose_matches)

    if total_units < 2:
        return False, f"Insufficient monosaccharide units, found {total_units}"

    # Define glycosidic linkage pattern (ether linkage between rings)
    glycosidic_pattern = Chem.MolFromSmarts("[C;!R]-O-[C;!R]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    if len(glycosidic_matches) < total_units - 1:
        return False, f"Insufficient glycosidic linkages for an oligosaccharide, found {len(glycosidic_matches)}"

    return True, "Molecule has multiple monosaccharide units connected by glycosidic linkages"