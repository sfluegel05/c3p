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

    # Normalize molecule (protonation states, tautomers)
    # This can be complex; for simplicity, let's assume input SMILES are standardized

    # Define SMARTS patterns for each B vitamin
    patterns = []

    # Vitamin B1 (Thiamine)
    # Thiazolium ring connected to pyrimidine ring via methylene bridge
    b1_pattern = Chem.MolFromSmarts("C[n+;H]1csc(C)c1C-Cc2ncnc(N)[nH]2")
    patterns.append(('Vitamin B1 (Thiamine)', b1_pattern))

    # Vitamin B2 (Riboflavin)
    # Isoalloxazine ring system with ribitol side chain
    b2_pattern = Chem.MolFromSmarts("C[C@@H](O)[C@@H](O)[C@@H](O)COc1nc2c3c(n1)nc(=O)[nH]c3=Nc2c1cc(C)c(C)cc1N2")
    patterns.append(('Vitamin B2 (Riboflavin)', b2_pattern))

    # Vitamin B3 (Niacin and Niacinamide)
    # Nicotinic acid or nicotinamide
    b3_pattern1 = Chem.MolFromSmarts("n1ccccc1C(=O)[O,N]")
    b3_pattern2 = Chem.MolFromSmarts("n1ccccc1C(=O)N")
    patterns.append(('Vitamin B3 (Niacin)', b3_pattern1))
    patterns.append(('Vitamin B3 (Niacinamide)', b3_pattern2))

    # Vitamin B5 (Pantothenic acid)
    # Pantoic acid moiety linked to beta-alanine via amide bond
    b5_pattern = Chem.MolFromSmarts("CC(C)(O)C(C(=O)NCCC(=O)[O,N])O")
    patterns.append(('Vitamin B5 (Pantothenic acid)', b5_pattern))

    # Vitamin B6 (Pyridoxine, Pyridoxal, Pyridoxamine)
    # Substituted pyridine rings, including phosphorylated forms
    # Pyridoxine
    b6_pyridoxine_pattern = Chem.MolFromSmarts("COc1cc(CO)c(C)c(O)n1")
    # Pyridoxal
    b6_pyridoxal_pattern = Chem.MolFromSmarts("O=CC1=NC=C(CO)C(C)C1O")
    # Pyridoxamine
    b6_pyridoxamine_pattern = Chem.MolFromSmarts("NCc1cc(CO)c(C)c(O)n1")
    # Phosphorylated forms
    b6_phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)OCC1=NC=C(CO)C(C)C1[O,N]")
    patterns.append(('Vitamin B6 (Pyridoxine)', b6_pyridoxine_pattern))
    patterns.append(('Vitamin B6 (Pyridoxal)', b6_pyridoxal_pattern))
    patterns.append(('Vitamin B6 (Pyridoxamine)', b6_pyridoxamine_pattern))
    patterns.append(('Vitamin B6 Phosphate', b6_phosphate_pattern))

    # Vitamin B7 (Biotin)
    # Heterocyclic fused ring system with valeric acid side chain
    b7_pattern = Chem.MolFromSmarts("O=C1NC(=O)N2C[C@H](SC1)[C@H]2CCCCCC(=O)[O,N]")
    patterns.append(('Vitamin B7 (Biotin)', b7_pattern))

    # Vitamin B9 (Folic acid and derivatives)
    # Pteridine ring connected to p-aminobenzoic acid and glutamic acid
    b9_pattern = Chem.MolFromSmarts("Nc1nc2[nH]ccc2n1-CNc1ccc(cc1)C(=O)N[C@@H](CCC(=O)[O,N])C(=O)[O,N]")
    patterns.append(('Vitamin B9 (Folic acid)', b9_pattern))

    # Vitamin B12 (Cobalamin)
    # Corrin ring with central cobalt atom and dimethylbenzimidazole
    b12_pattern = Chem.MolFromSmarts("[Co].[C,c]1=C[C,c]2[C,c]=[C,c][C,c]=[C,c][C,c]=[C,c][C,c]=1[C,c]=[C,c]2")
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