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
    
    # Detect sugar rings (5 or 6 membered rings with oxygen)
    sugar_ring_pattern = Chem.MolFromSmarts("[C,O]1[C,O][C,O][C,O][C,O][C,O]1")
    sugar_rings = mol.GetSubstructMatches(sugar_ring_pattern)
    if len(sugar_rings) < 2:
        return False, "Insufficient sugar rings found"
    
    # Detect glycosidic bonds (ether linkages connecting sugar rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C,O;R]-O-[C,O;R]")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_bonds) < 1:
        return False, "No glycosidic bonds found"
    
    # Check for uronic acid units (sugar ring with carboxylic acid group)
    uronic_acid_pattern = Chem.MolFromSmarts("[C;R][C;R][C;R][C;R][C;R][O;R]-C(=O)[O;H1,-1]")
    uronic_acid_units = mol.GetSubstructMatches(uronic_acid_pattern)
    if len(uronic_acid_units) < 1:
        return False, "No uronic acid units found"
    
    # Check for glycosamine units (sugar ring with amino group)
    glycosamine_pattern = Chem.MolFromSmarts("[C;R][C;R][C;R][C;R][C;R][O;R]-[CH2][N]")
    glycosamine_units = mol.GetSubstructMatches(glycosamine_pattern)
    if len(glycosamine_units) < 1:
        return False, "No glycosamine units found"
    
    # Check for sulfate ester groups (-O-S(=O)(=O)-O-)
    sulfate_ester_pattern = Chem.MolFromSmarts("O-S(=O)(=O)-O")
    sulfate_groups = mol.GetSubstructMatches(sulfate_ester_pattern)
    if len(sulfate_groups) < 1:
        return False, "No sulfate ester groups found"
    
    # Assuming that the presence of these units indicates a mucopolysaccharide
    return True, "Contains polysaccharide chain with uronic acid, glycosamine units, and sulfate groups"

__metadata__ = {
    'chemical_class': {
        'name': 'mucopolysaccharide',
        'definition': 'Any of the group of polysaccharides composed of alternating units from uronic acids and glycosamines, and commonly partially esterified with sulfuric acid.',
    }
}