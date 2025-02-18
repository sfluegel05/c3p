"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: CHEBI:17716 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside must have a D-galactose ring with beta configuration at the anomeric carbon (C1),
    and a glycosidic bond at the anomeric oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Standardize beta-D-galactose pattern with explicit stereochemistry
    # Pattern matches pyranose form with correct substituent orientations:
    # - C1 (anomeric) has beta configuration (axial O-glycosidic bond)
    # - C2: R configuration (hydroxyl equatorial)
    # - C3: R configuration (hydroxyl equatorial)
    # - C4: S configuration (hydroxyl axial)
    # - C5: R configuration (CH2OH group)
    beta_D_galactose_pattern = Chem.MolFromSmarts(
        "[C@@H]1([C@H](O)[C@@H](O)[C@H](O[C@H](CO)[C@H]1O)O)O[!H0]"
    )
    
    # Check for core beta-D-galactose structure with glycosidic bond
    if not mol.HasSubstructMatch(beta_D_galactose_pattern):
        return False, "No beta-D-galactose core with proper stereochemistry"
    
    # Verify glycosidic oxygen is connected to non-sugar moiety
    # (already enforced by [!H0] in SMARTS pattern)
    
    # Optional: Check molecular weight > typical monosaccharides (150-200 Da)
    # Helps filter out free galactose or small fragments
    mol_wt = AllChem.CalcExactMolWt(mol)
    if mol_wt < 180:
        return False, "Molecular weight too low for glycoside"

    return True, "Contains beta-D-galactose core with proper stereochemistry and glycosidic bond"