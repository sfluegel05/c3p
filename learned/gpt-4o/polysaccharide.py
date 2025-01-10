"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide contains more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")
    
    # Define SMARTS pattern for pyranose and furanose structures, common in polysaccharides
    pyranose = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H]1")
    furanose = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@@H]1")
    
    # Counting monosaccharide units
    pyranose_matches = len(mol.GetSubstructMatches(pyranose))
    furanose_matches = len(mol.GetSubstructMatches(furanose))
    total_monosaccharides = pyranose_matches + furanose_matches

    if total_monosaccharides <= 10:
        return (False, f"Only found {total_monosaccharides} monosaccharide units, require more than 10")

    # Check for glycosidic linkages
    glycosidic_linkage = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H]1")
    if not mol.HasSubstructMatch(glycosidic_linkage):
        return (False, "No glycosidic linkages found")

    # Confirm sufficient numbers and linkages
    return (True, "Valid polysaccharide with necessary glycosidic linkages and sufficient monosaccharide units")