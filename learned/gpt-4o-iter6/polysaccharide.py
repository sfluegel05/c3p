"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide consists of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for hexoses - commonly seen monosaccharides
    # C1(OC)OC(C)OCC1 is a common simplification for hexopyranoses
    hexose_pattern = Chem.MolFromSmarts("C1(O)[C@H](O)[C@H](O)[C@@H](O[C@H]1)CO") 
    # Acetylated or amine groups in common oligosaccharide repeating units (e.g., acetylglucosamine)
    derivative_pattern = Chem.MolFromSmarts("NC(=O)C")

    # Match and count hexose-derivatives
    hexose_count = len(mol.GetSubstructMatches(hexose_pattern))
    derivative_count = len(mol.GetSubstructMatches(derivative_pattern))

    # Sum the relevant sugar-like components
    sugar_unit_count = hexose_count + derivative_count

    if sugar_unit_count <= 10:
        return False, f"Only found {sugar_unit_count} monosaccharide-like units, requiring more than 10 for classification as a polysaccharide."

    # Define pattern for glycosidic linkages
    glyco_linkage_pattern = Chem.MolFromSmarts("[C@H]([C@@H]O[C@H])")  # general linkage presence
    glyco_link_count = len(mol.GetSubstructMatches(glyco_linkage_pattern))

    # Verify sufficient glycosidic bonds
    if glyco_link_count < sugar_unit_count - 1:
        return False, f"Detected {glyco_link_count} glycosidic bonds; should be at least {sugar_unit_count - 1} to connect the units."

    # Criteria satisfied
    return True, "Contains adequate monosaccharide units and glycosidic linkages to be classified as a polysaccharide."