"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    The B vitamins include B1 (thiamine), B2 (riboflavin), B3 (niacin), B5 (pantothenic acid),
    B6 (pyridoxine, pyridoxal, pyridoxamine), B7 (biotin), B9 (folic acid), and B12 (cobalamin).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for each B vitamin
    patterns = []

    # Vitamin B1 (thiamine)
    # Thiazole ring connected to pyrimidine ring via methylene bridge
    b1_pattern = Chem.MolFromSmarts("Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1")
    patterns.append(('Vitamin B1 (Thiamine)', b1_pattern))

    # Vitamin B2 (riboflavin)
    # Isoalloxazine ring system
    b2_pattern = Chem.MolFromSmarts("Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H](O)[C@H](O)[C@H](O)CO)c2cc1C")
    patterns.append(('Vitamin B2 (Riboflavin)', b2_pattern))

    # Vitamin B3 (niacin)
    # Pyridine ring with carboxylic acid at position 3
    b3_pattern = Chem.MolFromSmarts("c1cc(cnc1)C(=O)[O,H]")
    patterns.append(('Vitamin B3 (Niacin)', b3_pattern))

    # Vitamin B5 (pantothenic acid)
    # Pantoic acid moiety linked to beta-alanine via amide bond
    b5_pattern = Chem.MolFromSmarts("CC(C)(O)C(C(=O)NCCC(=O)[O,H])O")
    patterns.append(('Vitamin B5 (Pantothenic acid)', b5_pattern))

    # Vitamin B6 (pyridoxine, pyridoxal, pyridoxamine)
    # Substituted pyridine rings
    # Pyridoxine
    b6_pyridoxine_pattern = Chem.MolFromSmarts("COc1cc(CO)c(C)c(O)n1")
    patterns.append(('Vitamin B6 (Pyridoxine)', b6_pyridoxine_pattern))
    # Pyridoxal
    b6_pyridoxal_pattern = Chem.MolFromSmarts("O=CC1=NC=C(CO)C(C)C1O")
    patterns.append(('Vitamin B6 (Pyridoxal)', b6_pyridoxal_pattern))
    # Pyridoxamine
    b6_pyridoxamine_pattern = Chem.MolFromSmarts("NCc1cc(CO)c(C)c(O)n1")
    patterns.append(('Vitamin B6 (Pyridoxamine)', b6_pyridoxamine_pattern))

    # Vitamin B7 (biotin)
    # Fused ureido and tetrahydrothiophene rings
    b7_pattern = Chem.MolFromSmarts("O=C1NC(=O)N2[C@@H](CS[C@@H]2CCCCC(=O)[O,H])[C@H]1")
    patterns.append(('Vitamin B7 (Biotin)', b7_pattern))

    # Vitamin B9 (folic acid)
    # Pteridine ring system connected to p-aminobenzoic acid and glutamic acid
    b9_pattern = Chem.MolFromSmarts("Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)[O,H])C(=O)[O,H])Nc2c(=O)[nH]1")
    patterns.append(('Vitamin B9 (Folic acid)', b9_pattern))

    # Vitamin B12 (cobalamin)
    # Complex corrin ring with cobalt ion
    # Simplified pattern to check for cobalt atom in corrin ring
    b12_pattern = Chem.MolFromSmarts("[Co]")  # Checks for presence of cobalt
    patterns.append(('Vitamin B12 (Cobalamin)', b12_pattern))

    # Check for matches
    for name, pattern in patterns:
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches {name} pattern"

    return False, "Does not match any B vitamin patterns"

__metadata__ = {
    'chemical_class': {
        'name': 'B vitamin',
        'definition': 'Any member of the group of eight water-soluble vitamins originally thought to be a single compound (vitamin B) that play important roles in cell metabolism. The group comprises vitamin B1, B2, B3, B5, B6, B7, B9, and B12.',
    }
}