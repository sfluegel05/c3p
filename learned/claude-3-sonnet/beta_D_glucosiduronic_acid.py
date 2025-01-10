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

    # SMARTS patterns for beta-D-glucuronic acid core
    # Pattern matches the full beta-D-glucuronic acid core with correct stereochemistry
    # [C@@H]1 at anomeric carbon (beta), followed by the ring with correct D-configuration
    # and carboxylic acid at C6
    glucuronic_core = Chem.MolFromSmarts(
        '[OD2][C@@H]1[C@H]([OH1])[C@@H]([OH1])[C@H]([OH1])[C@@H](C(=O)[OH1])O1'
    )
    
    # Alternative pattern that matches the core when it's part of a glycoside
    glucuronic_glycoside = Chem.MolFromSmarts(
        '[OD2][C@@H]1[C@H]([OH1])[C@@H]([OH1])[C@H]([OH1])[C@@H](C(=O)[OH1])O1'
    )

    # Check for carboxylic acid at C6
    carboxyl = Chem.MolFromSmarts('C(=O)[OH1]')
    if not mol.HasSubstructMatch(carboxyl):
        return False, "No carboxylic acid group found"

    # Check for pyranose ring with correct stereochemistry
    if not (mol.HasSubstructMatch(glucuronic_core) or mol.HasSubstructMatch(glucuronic_glycoside)):
        return False, "No beta-D-glucuronic acid core found with correct stereochemistry"

    # Check for proper number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Insufficient oxygen atoms for glucuronic acid structure"

    # Check for glycosidic linkage
    # The anomeric carbon should be connected to two oxygens
    anomeric_pattern = Chem.MolFromSmarts('[OD2][C@@H]([OH0])[C@H]([OH1])')
    if not mol.HasSubstructMatch(anomeric_pattern):
        return False, "No proper glycosidic linkage found"

    # Verify hydroxyl groups at correct positions (excluding the carboxylic acid)
    hydroxyl_pattern = Chem.MolFromSmarts('[C@@H]([OH1])[C@H]([OH1])[C@@H]([OH1])')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing required hydroxyl groups"

    # Additional check for six-membered ring
    ring_info = mol.GetRingInfo()
    if not any(size == 6 for size in ring_info.RingSizes()):
        return False, "No six-membered ring found"

    return True, "Contains beta-D-glucuronic acid core with glycosidic linkage"