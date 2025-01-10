"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide has repeated units of monosaccharides linked by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycosidic linkage pattern
    glycosidic_pattern = Chem.MolFromSmarts("O[C@@H]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"
    
    # Look for repeated cyclic monosaccharide units
    cyclic_monosaccharide_pattern = Chem.MolFromSmarts("C1([C@H](O)C([C@H](O)[C@H](O)[C@@H]1O)O)")
    repeat_units = mol.GetSubstructMatches(cyclic_monosaccharide_pattern)
    if len(repeat_units) < 2:
        return False, "Too few repeating units for polysaccharide"
    
    # Check number of hydroxyl groups
    num_oh = sum(atom.GetSymbol() == 'O' and len(atom.GetNeighbors()) == 1 for atom in mol.GetAtoms())
    if num_oh < 5:
        return False, "Too few hydroxyl groups"

    # Check molecular weight - polysaccharides are often large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for polysaccharide"

    return True, "Contains multiple repeated monosaccharide units linked by glycosidic bonds"