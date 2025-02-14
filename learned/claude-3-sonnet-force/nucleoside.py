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
    sugar_pattern = Chem.MolFromSmarts("[OX2][CX4]([CX4][CX4][CX4][CH2][OX2])")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No ribose or deoxyribose sugar found"

    # Look for nucleobase
    nucleobase_patterns = [
        Chem.MolFromSmarts("[*r3,r4,r5,r6]~[*r3,r4,r5,r6]~[*r3,r4,r5,r6]~[*r3,r4,r5,r6]~[NX3]"),  # Fused ring system with nitrogen
        Chem.MolFromSmarts("[NX3][CX3]([NX3])[NX3]"),  # Cytosine-like
        Chem.MolFromSmarts("[NX3][CX3]([NX3])[OX2]"),  # Uracil-like
        # Add more patterns for other nucleobases
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns):
        return False, "No nucleobase found"

    # Check for N-glycosidic bond between sugar and nucleobase
    for sugar_match in sugar_matches:
        sugar_atom = mol.GetAtomWithIdx(sugar_match[-1])
        for neighbor_atom in sugar_atom.GetNeighbors():
            if neighbor_atom.GetAtomicNum() == 7 and neighbor_atom.IsInRing():  # Nitrogen atom in a ring (nucleobase)
                return True, "Contains a nucleobase N-glycosidically linked to a ribose or deoxyribose sugar"

    return False, "No N-glycosidic bond found between sugar and nucleobase"