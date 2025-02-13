"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: mucopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating units of uronic acids and glycosamines,
    and commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycosidic bonds (O-linked sugars)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C;!$(C=O)]1O[C;!$(C=O)][C;!$(C=O)][C;!$(C=O)][C;!$(C=O)][O][C;!$(C=O)]1")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_bonds:
        return False, "No glycosidic bonds found"
    
    # Check for uronic acid units (sugar ring with carboxylic acid)
    uronic_acid_pattern = Chem.MolFromSmarts("[C@@H]1([O])[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1C(=O)[O;H1,-]")
    uronic_acid_units = mol.GetSubstructMatches(uronic_acid_pattern)
    if not uronic_acid_units:
        return False, "No uronic acid units found"
    
    # Check for glycosamine units (sugar ring with amino group)
    glycosamine_pattern = Chem.MolFromSmarts("[C@@H]1([O])[C@@H](O)[C@H](O)[C@@H](O)[C@H](N)[C@@H]1O")
    glycosamine_units = mol.GetSubstructMatches(glycosamine_pattern)
    if not glycosamine_units:
        return False, "No glycosamine units found"
    
    # Check for sulfate ester groups (-O-S(=O)(=O)-O-)
    sulfate_ester_pattern = Chem.MolFromSmarts("O[S](=O)(=O)[O]")
    sulfate_groups = mol.GetSubstructMatches(sulfate_ester_pattern)
    if not sulfate_groups:
        return False, "No sulfate ester groups found"
    
    # Check for alternating units (simplified check)
    if len(uronic_acid_units) < 1 or len(glycosamine_units) < 1:
        return False, "Does not contain both uronic acid and glycosamine units"
    
    # Assuming that the presence of these units indicates alternating pattern
    return True, "Contains polysaccharide chain with uronic acid and glycosamine units and sulfate groups"

__metadata__ = {
    'chemical_class': {
        'name': 'mucopolysaccharide',
        'definition': 'Any of the group of polysaccharides composed of alternating units from uronic acids and glycosamines, and commonly partially esterified with sulfuric acid.',
    }
}