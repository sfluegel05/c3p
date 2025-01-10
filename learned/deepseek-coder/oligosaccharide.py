"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound in which monosaccharide units are joined by glycosidic linkages.

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

    # Look for glycosidic linkage pattern (more specific pattern)
    glycosidic_pattern = Chem.MolFromSmarts("[C;H1,H2][O][C;H1,H2]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Look for multiple hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 4:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 4"

    # Check for multiple sugar units (more flexible pattern)
    sugar_unit_pattern = Chem.MolFromSmarts("[C;H1,H2][C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    sugar_unit_matches = mol.GetSubstructMatches(sugar_unit_pattern)
    if len(sugar_unit_matches) < 2:
        return False, f"Found {len(sugar_unit_matches)} sugar units, need at least 2"

    # Check molecular weight - oligosaccharides typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for oligosaccharide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:
        return False, "Too few carbons for oligosaccharide"
    if o_count < 6:
        return False, "Too few oxygens for oligosaccharide"

    # Check for common sugar configurations
    sugar_config_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@H](O)[C@H](O)")
    if not mol.HasSubstructMatch(sugar_config_pattern):
        return False, "No typical sugar configuration found"

    return True, "Contains multiple sugar units joined by glycosidic linkages"