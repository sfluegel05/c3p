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
    The B vitamins include B1 (thiamine), B2 (riboflavin), B3 (niacin),
    B5 (pantothenic acid), B6 (pyridoxine, pyridoxal, pyridoxamine),
    B7 (biotin), B9 (folic acid), and B12 (cobalamin).

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

    # Standardize molecule (e.g., remove salts, normalize)
    mol = Chem.RemoveHs(mol)

    # Define substructure patterns for each B vitamin

    patterns = []

    # Vitamin B1 (Thiamine) - Thiazolium ring connected to pyrimidine ring
    b1_pattern = Chem.MolFromSmarts("c1ncc[nH]1CCc1[n+](C)sc(C)c1")
    patterns.append(('Vitamin B1 (Thiamine)', b1_pattern))

    # Vitamin B2 (Riboflavin) - Isoalloxazine ring with ribitol side chain
    b2_pattern = Chem.MolFromSmarts("C1=CC2=C3C(=O)N=CNC3=NC(=O)N2C(=C1)N[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O")
    patterns.append(('Vitamin B2 (Riboflavin)', b2_pattern))

    # Vitamin B3 (Niacin and Niacinamide) - Nicotinic acid and nicotinamide
    b3_pattern = Chem.MolFromSmarts("n1ccccc1C(=O)[O,N]")
    patterns.append(('Vitamin B3 (Niacin/Niacinamide)', b3_pattern))

    # Vitamin B5 (Pantothenic acid)
    b5_pattern = Chem.MolFromSmarts("CC(C)(O)C(O)C(=O)NCCC(=O)O")
    patterns.append(('Vitamin B5 (Pantothenic acid)', b5_pattern))

    # Vitamin B6 (Pyridoxine, Pyridoxal, Pyridoxamine)
    b6_pattern = Chem.MolFromSmarts("OCc1ccn(C)[nH]c1C")
    patterns.append(('Vitamin B6 (Pyridoxine family)', b6_pattern))

    # Vitamin B7 (Biotin) - Imidazolidone ring fused to a tetrahydrothiophene ring
    b7_pattern = Chem.MolFromSmarts("O=C1NC(=O)[C@@H]2CSC[C@H]2N1")
    patterns.append(('Vitamin B7 (Biotin)', b7_pattern))

    # Vitamin B9 (Folic acid and derivatives)
    b9_pattern = Chem.MolFromSmarts("Nc1nc2ncnc(N)c2n1-c1ccc(cc1)C(=O)N[C@H](CCC(=O)[O,N])C(=O)[O,N]")
    patterns.append(('Vitamin B9 (Folic acid)', b9_pattern))

    # Vitamin B12 (Cobalamin) - Corrin ring with cobalt
    b12_pattern = Chem.MolFromSmarts("[Co]")
    patterns.append(('Vitamin B12 (Cobalamin)', b12_pattern))

    # Check for matches
    for name, pattern in patterns:
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches {name}"

    # Additional check: Molecules containing cobalt for B12
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 27:
            return True, "Contains cobalt atom indicative of Vitamin B12"

    # Check for phosphorylated forms of B vitamins
    phosphorylated = False
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])[O-]")
    if mol.HasSubstructMatch(phosphate_pattern):
        phosphorylated = True

    # Specific check for phosphorylated B6 vitamins
    b6_phosphate_pattern = Chem.MolFromSmarts("O[P](=O)([O-])[O-]COc1ccn(C)[nH]c1C")
    if mol.HasSubstructMatch(b6_phosphate_pattern):
        return True, "Matches Vitamin B6 phosphorylated form"

    # Specific check for Thiamine pyrophosphate (B1 active form)
    b1_pyrophosphate_pattern = Chem.MolFromSmarts("C[n+](C)c1ncccn1CCOP(=O)([O-])OP(=O)([O-])[O-]")
    if mol.HasSubstructMatch(b1_pyrophosphate_pattern):
        return True, "Matches Vitamin B1 pyrophosphate form"

    # Fallback similarity check with known B vitamins
    known_b_vitamins_smiles = [
        # SMILES of known B vitamins
        "C[n+](C)c1ncccn1CCc1csc(C)[n+]1C",  # Thiamine
        "Cc1c2nc3nc(c(=O)[nH]c3=O)n(C)c2nc1N[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",  # Riboflavin
        "OC(=O)c1cccnc1",  # Niacin
        "OC(=O)CNC(CC(O)C(C)(C)O)C(=O)O",  # Pantothenic acid
        "OCc1cc(C)n(C)c(=O)c1O",  # Pyridoxal
        "C[C@@H](O)[C@H](O)[C@H](COP(=O)([O-])[O-])O",  # Phosphorylated ribose
        "O=C1NC(=O)[C@H]2CSC[C@H]2N1",  # Biotin
        "Nc1nc2ncnc(N)c2n1-c1ccc(cc1)C(=O)O",  # Folic acid
        # ... Add more as needed
    ]
    known_b_vitamins = [Chem.MolFromSmiles(smi) for smi in known_b_vitamins_smiles if Chem.MolFromSmiles(smi) is not None]

    fp_query = AllChem.GetMorganFingerprint(mol, 2)

    for idx, vitamin_mol in enumerate(known_b_vitamins):
        fp_vitamin = AllChem.GetMorganFingerprint(vitamin_mol, 2)
        similarity = Chem.DataStructs.TanimotoSimilarity(fp_query, fp_vitamin)
        if similarity > 0.7:
            return True, f"Similar to known B vitamin with Tanimoto similarity {similarity:.2f}"

    return False, "Does not match any B vitamin patterns"

__metadata__ = {
    'chemical_class': {
        'name': 'B vitamin',
        'definition': 'Any member of the group of eight water-soluble vitamins originally thought to be a single compound (vitamin B) that play important roles in cell metabolism. The group comprises vitamin B1, B2, B3, B5, B6, B7, B9, and B12.',
    }
}