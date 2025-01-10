"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: CHEBI:16761 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleoside with one or more phosphate groups attached to the sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sugar pattern (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]([C@H](O)[C@@H]1O)")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar pattern found"

    # Look for nucleobase pattern (more general pattern)
    nucleobase_pattern = Chem.MolFromSmarts("[nX3]1[cX3][cX3][cX3][cX3]1")
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase pattern found"

    # Look for phosphate groups (more general pattern)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"

    # Check molecular weight - nucleoside phosphates typically >150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for nucleoside phosphate"

    # Count phosphorous atoms
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 1:
        return False, "Must have at least one phosphorous atom"

    return True, "Contains sugar, nucleobase, and at least one phosphate group"