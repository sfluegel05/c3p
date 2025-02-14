"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
"""

from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    The B vitamins include vitamin B1 (thiamine), B2 (riboflavin), B3 (niacin),
    B5 (pantothenic acid), B6 (pyridoxine, pyridoxal, pyridoxamine),
    B7 (biotin), B9 (folic acid), and B12 (cobalamin).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a B vitamin, False otherwise
        str: Reason for classification
    """

    # Convert input SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for each B vitamin
    b_vitamin_patterns = [
        # Vitamin B1 (Thiamine) - thiazolium ring connected to a pyrimidine ring
        (Chem.MolFromSmarts('c1[n+](ccs1)CCc1ncccn1'), 'Vitamin B1 (Thiamine)'),
        # Vitamin B2 (Riboflavin) - isoalloxazine ring system
        (Chem.MolFromSmarts('c1cc2nc3c([nH]c(=O)[nH]c3nc2c(c1)N)[C@H]1O[C@H](CO)[C@H](O)[C@H]1O'), 'Vitamin B2 (Riboflavin)'),
        # Vitamin B3 (Niacin) - nicotinic acid and nicotinamide
        (Chem.MolFromSmarts('c1cccnc1C(=O)[O,NH2]'), 'Vitamin B3 (Niacin)'),
        # Vitamin B5 (Pantothenic acid) - beta-alanine linked to pantoic acid
        (Chem.MolFromSmarts('CC(C)(CO)C(=O)NCCC(=O)O'), 'Vitamin B5 (Pantothenic acid)'),
        # Vitamin B6 (Pyridoxine, Pyridoxal, Pyridoxamine)
        (Chem.MolFromSmarts('c1cc(CO)c(C)c(O)c1'), 'Vitamin B6 (Pyridoxine/Pyridoxal/Pyridoxamine)'),
        # Vitamin B7 (Biotin) - ureido ring fused with a tetrahydrothiophene ring
        (Chem.MolFromSmarts('O=C1NC(=O)N2C[C@@H](SC1)[C@]2([H])'), 'Vitamin B7 (Biotin)'),
        # Vitamin B9 (Folic acid and derivatives) - pteridine ring linked to p-aminobenzoic acid and glutamic acid
        (Chem.MolFromSmarts('Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CC(O)=O)C(=O)O)nc2c(=O)[nH]1'), 'Vitamin B9 (Folic acid)'),
        # Vitamin B12 (Cobalamin) - corrin ring with central cobalt atom
        (Chem.MolFromSmarts('[Cobalt]'), 'Vitamin B12 (Cobalamin)'),
    ]

    # Check if the molecule matches any B vitamin pattern
    for pattern, vit_name in b_vitamin_patterns:
        if pattern is None:
            continue  # Skip if the pattern is invalid
        if mol.HasSubstructMatch(pattern):
            return True, f"Molecule matches {vit_name}"

    # Special check for cobalamin derivatives (Vitamin B12)
    if any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms()):
        return True, "Molecule contains cobalt atom characteristic of Vitamin B12 (Cobalamin)"

    return False, "Molecule does not match any known B vitamin"