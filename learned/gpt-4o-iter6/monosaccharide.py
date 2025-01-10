"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon and oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Monosaccharides must have at least 3 carbons
    if c_count < 3:
        return False, "Too few carbon atoms"

    # Check for hydroxyl groups (at least two -OH groups)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Less than 2 hydroxyl groups"

    # Check for carbonyl group (aldehyde C=O or ketone C=O)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group (aldehyde or ketone)"

    # Check that it's a single unit (no glycosidic bond)
    # Glycosidic linkage typically involves anomeric carbon connecting via an ether linkage
    anomeric_ether_pattern = Chem.MolFromSmarts("[C@H](O)([CH0,CH1,CH2,CH3,O])O")
    if mol.HasSubstructMatch(anomeric_ether_pattern):
        return False, "Contains glycosidic connection"

    # Checking for deoxy sugars, which still have carbonyl group
    deoxy_pattern = Chem.MolFromSmarts("[CH]C(=O)")
    if mol.HasSubstructMatch(deoxy_pattern):
        if not mol.HasSubstructMatch(carbonyl_pattern):
            return False, "Lacks necessary structural features of monosaccharide"

    return True, "Matches structure of a monosaccharide with polyhydroxyl groups and a carbonyl group"