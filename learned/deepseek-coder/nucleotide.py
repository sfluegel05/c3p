"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for nucleoside (sugar + base) and phosphate
    sugar_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]([C@H](O)[C@@H]1O)")  # Ribose or deoxyribose sugar
    base_pattern = Chem.MolFromSmarts("[nH]1cnc2c1nc[nH]2")  # Adenine, guanine, cytosine, thymine, or uracil
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")  # Phosphate group

    # Check for nucleoside (sugar + base)
    if not mol.HasSubstructMatch(sugar_pattern) or not mol.HasSubstructMatch(base_pattern):
        return False, "No nucleoside (sugar + base) found"

    # Check for phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check if phosphate is attached to sugar (3' or 5' position)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    phosphate_attached = False

    for phosphate_match in phosphate_matches:
        for sugar_match in sugar_matches:
            # Check if phosphate is attached to sugar
            if any(atom_idx in sugar_match for atom_idx in phosphate_match):
                phosphate_attached = True
                break
        if phosphate_attached:
            break

    if not phosphate_attached:
        return False, "Phosphate group not attached to sugar"

    return True, "Contains nucleoside with phosphate group attached to sugar"