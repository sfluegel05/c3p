"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for monosaccharide units (pyranose and furanose rings)
    # Pyranose: six-membered ring with 5 carbons and 1 oxygen
    pyranose = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)O1")
    # Furanose: five-membered ring with 4 carbons and 1 oxygen
    furanose = Chem.MolFromSmarts("C1C(O)C(O)C(O)O1")

    # Find all monosaccharide units
    pyranose_matches = mol.GetSubstructMatches(pyranose)
    furanose_matches = mol.GetSubstructMatches(furanose)

    total_units = len(pyranose_matches) + len(furanose_matches)
    if total_units !=4:
        return False, f"Found {total_units} monosaccharide units, expected 4"

    # Optional: Check for glycosidic linkages
    # Glycosidic linkage pattern: anomeric carbon connected via oxygen to another sugar
    glycosidic_linkage = Chem.MolFromSmarts("[C;R][O][C;R]")
    linkage_matches = mol.GetSubstructMatches(glycosidic_linkage)
    if len(linkage_matches) < 3:
        return False, "Not enough glycosidic linkages between monosaccharide units"

    # Check molecular weight range for typical tetrasaccharides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.2f} not typical for a tetrasaccharide"

    return True, "Molecule contains 4 monosaccharide units connected via glycosidic linkages"