"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:36314 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide consists of a trisaccharide repeating unit (two heptose units and octulosonic acid),
    oligosaccharide side chains, and 3-hydroxytetradecanoic acid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for trisaccharide repeating unit
    trisaccharide_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@@H]([C@H]([C@@H]1O)O)O)O[C@H]2[C@@H]([C@H]([C@H]([C@H]2O)O[C@H]3[C@@H]([C@H]([C@H]([C@@H](O3)O)O)O)O")
    if not mol.HasSubstructMatch(trisaccharide_pattern):
        return False, "Trisaccharide repeating unit not found"

    # Look for oligosaccharide side chains
    oligosaccharide_pattern = Chem.MolFromSmarts("[OX2][C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O[C@H]2[C@@H]([C@H]([C@H]([C@@H](O2)O)O)O")
    if not mol.HasSubstructMatch(oligosaccharide_pattern):
        return False, "Oligosaccharide side chains not found"

    # Look for 3-hydroxytetradecanoic acid units
    hydroxy_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(O)CC(O)=O")
    if not mol.HasSubstructMatch(hydroxy_acid_pattern):
        return False, "3-hydroxytetradecanoic acid units not found"

    # Check molecular weight - lipopolysaccharides typically >2000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 2000:
        return False, "Molecular weight too low for lipopolysaccharide"

    return True, "Contains trisaccharide repeating unit, oligosaccharide side chains, and 3-hydroxytetradecanoic acid units"