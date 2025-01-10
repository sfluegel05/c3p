"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is formed by the condensation of any substance with 
    beta-D-glucuronic acid to form a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for beta-D-glucuronic acid core with correct stereochemistry
    # Matches the complete core structure including the carboxylic acid
    glucuronic_pattern = Chem.MolFromSmarts(
        '[OD2][C@@H]1[C@H]([OH1])[C@@H]([OH1])[C@H]([OH1])[C@@H](C(=O)[OH1])O1'
    )
    
    # Pattern for when the core is part of a glycoside (O-substituted)
    glucuronic_glycoside = Chem.MolFromSmarts(
        '[OD2][C@@H]1[C@H]([OH1])[C@@H]([OH1])[C@H]([OH1])[C@@H](C(=O)[OH1])O1'
    )

    # Check for the presence of carboxylic acid group
    carboxyl = Chem.MolFromSmarts('C(=O)[OH1]')
    if not mol.HasSubstructMatch(carboxyl):
        return False, "No carboxylic acid group found"

    # Check for the presence of pyranose ring with correct stereochemistry
    has_core = mol.HasSubstructMatch(glucuronic_pattern)
    has_glycoside = mol.HasSubstructMatch(glucuronic_glycoside)
    
    if not (has_core or has_glycoside):
        return False, "No beta-D-glucuronic acid core found with correct stereochemistry"

    # Check for proper glycosidic linkage at anomeric carbon
    anomeric_o = Chem.MolFromSmarts('[OD2][C@@H]1O[C@@H]')
    if not mol.HasSubstructMatch(anomeric_o):
        return False, "No proper glycosidic linkage found at anomeric carbon"

    # Verify presence of required hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts('[C@@H]([OH1])[C@H]([OH1])[C@@H]([OH1])')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing required hydroxyl groups"

    # Count oxygen atoms (should have at least 6 for glucuronic acid)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Insufficient oxygen atoms for glucuronic acid structure"

    # Check for 6-membered ring containing the glucuronic acid core
    rings = mol.GetRingInfo()
    if rings.NumRings() == 0:
        return False, "No rings found in molecule"

    # Additional check for correct substitution pattern
    pyranose_pattern = Chem.MolFromSmarts('[OD2][C@@H]1[C@H][C@@H][C@H][C@@H](C(=O)[OH1])O1')
    if not mol.HasSubstructMatch(pyranose_pattern):
        return False, "Incorrect substitution pattern for beta-D-glucuronic acid"

    return True, "Contains beta-D-glucuronic acid core with correct stereochemistry and glycosidic linkage"