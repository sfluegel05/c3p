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

    # Look for nucleoside pattern (sugar + nucleobase)
    nucleoside_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(nucleoside_pattern):
        return False, "No nucleoside pattern found (sugar + nucleobase)"

    # Look for phosphate groups (P with oxygens)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX2])([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"

    # Check molecular weight - nucleoside phosphates typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for nucleoside phosphate"

    # Count phosphorous atoms
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 1:
        return False, "Must have at least one phosphorous atom"

    return True, "Contains nucleoside pattern with at least one phosphate group"