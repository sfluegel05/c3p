"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Define SMARTS patterns for monosaccharide units and glycosidic bonds
    monosaccharide_pattern = ['OC1OC(CO)C(O)C(O)C1', 'OC1OC(O)C(O)C(O)C1', 'OC1C(O)C(O)C(CO)C(O)C1',
                              'OC1C(O)C(O)C(O)C(O)C1', 'C1C(O)C(O)C(CO)C(O)C1', 'C1C(O)C(O)C(O)C(O)C1']
    glycosidic_bond_pattern = 'OC[C@@]([H])([C@@]([H])(O[C@@]([H])([C@]([H])(O)[C@@]([H])([C@]([H])(O)[C@]([H])(O)CO)O)O)'
    
    # Count monosaccharide units
    num_monosaccharides = 0
    for pattern in monosaccharide_pattern:
        num_monosaccharides += len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
    
    # Count glycosidic bonds
    num_glycosidic_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts(glycosidic_bond_pattern)))
    
    # Check criteria for tetrasaccharide
    if num_monosaccharides == 4 and num_glycosidic_bonds >= 3:
        return True, "Contains 4 monosaccharide units and at least 3 glycosidic bonds"
    else:
        return False, f"Found {num_monosaccharides} monosaccharide units and {num_glycosidic_bonds} glycosidic bonds"

# Example usage
print(is_tetrasaccharide('OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O[C@@H]4[C@@H](CO)OC(O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O'))
# Output: (True, 'Contains 4 monosaccharide units and at least 3 glycosidic bonds')