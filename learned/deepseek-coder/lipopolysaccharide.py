"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide is a complex molecule consisting of a trisaccharide repeating unit
    (two heptose units and octulosonic acid) with oligosaccharide side chains and
    3-hydroxytetradecanoic acid units.

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

    # Look for trisaccharide repeating unit pattern (general pattern for three connected sugar units)
    trisaccharide_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O)[C@@H](O)[C@@H](O1)O)[C@H]2([C@H](O)[C@H](O)[C@@H](O)[C@@H](O2)O)[C@H]3([C@H](O)[C@H](O)[C@@H](O)[C@@H](O3)O)")
    if not mol.HasSubstructMatch(trisaccharide_pattern):
        return False, "No trisaccharide repeating unit found"

    # Look for oligosaccharide side chains (multiple sugar units)
    oligosaccharide_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O)[C@@H](O)[C@@H](O1)O)")
    oligosaccharide_matches = mol.GetSubstructMatches(oligosaccharide_pattern)
    if len(oligosaccharide_matches) < 3:
        return False, f"Found {len(oligosaccharide_matches)} oligosaccharide units, need at least 3"

    # Look for 3-hydroxytetradecanoic acid units (general pattern for long carbon chain with hydroxyl and carboxyl groups)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]([OH])[CX3](=[OX1])[OX2H]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Found {len(fatty_acid_matches)} 3-hydroxytetradecanoic acid units, need at least 1"

    # Check molecular weight - lipopolysaccharides typically >1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for lipopolysaccharide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for lipopolysaccharide"
    if o_count < 8:
        return False, "Too few oxygens for lipopolysaccharide"

    return True, "Contains trisaccharide repeating unit with oligosaccharide side chains and 3-hydroxytetradecanoic acid units"