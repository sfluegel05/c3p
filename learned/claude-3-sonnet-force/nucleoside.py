"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:18010 nucleoside
An N-glycosyl compound that has both a nucleobase, normally adenine, guanine, xanthine, thymine, cytosine or uracil, and either a ribose or deoxyribose as functional parents.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ribose or deoxyribose sugar
    sugar_patterns = [
        Chem.MolFromSmarts("[OX2]C[CX4]([OX2])[CX4]([CX4]([OX2])[CX4][OX2])[CH2][OX2]"),
        Chem.MolFromSmarts("[OX2]C[CX4]([OX2])[CX4]([CX4]([OX2])[CX4][OX2])[CH2][OX2H]")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns):
        return False, "No ribose or deoxyribose sugar found"

    # Look for nucleobase
    nucleobase_pattern = Chem.MolFromSmarts("[*r3,r4,r5,r6]~[*r3,r4,r5,r6]~[*r3,r4,r5,r6]~[*r3,r4,r5,r6]~[NX3]")
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"

    # Check for N-glycosidic bond between sugar and nucleobase
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][CX4][CX4][NX3]")
    if not AllChem.FindMolChiralCenters(mol, glycosidic_bond_pattern):
        return False, "No N-glycosidic bond found between sugar and nucleobase"

    # Additional checks
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 3:
        return False, "Too few nitrogen atoms for a nucleoside"

    return True, "Contains a nucleobase N-glycosidically linked to a ribose or deoxyribose sugar"