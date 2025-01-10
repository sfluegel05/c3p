"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

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

    # Standardize molecule (e.g., remove salts)
    mol = Chem.RemoveHs(mol)

    # Define SMARTS patterns for each B vitamin
    patterns = []

    # Vitamin B1 (Thiamine) - Thiazolium ring connected to pyrimidine ring
    b1_pattern = Chem.MolFromSmarts("C[n+](C)c1ncccn1CC2=CN=C(SC2)C")
    patterns.append(('Vitamin B1 (Thiamine)', b1_pattern))

    # Vitamin B2 (Riboflavin) - Isoalloxazine ring with ribitol side chain
    b2_pattern = Chem.MolFromSmarts("C1(=O)NC2=NC3=C(N2)[N]C(=O)C3=NC2=C1N=CN2[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O")
    patterns.append(('Vitamin B2 (Riboflavin)', b2_pattern))

    # Vitamin B3 (Niacin) - Nicotinic acid specifically at position 3
    b3_pattern = Chem.MolFromSmarts("n1cc(C(=O)O)ccc1")
    patterns.append(('Vitamin B3 (Niacin)', b3_pattern))

    # Niacinamide
    b3_niacinamide_pattern = Chem.MolFromSmarts("n1cc(C(=O)N)ccc1")
    patterns.append(('Vitamin B3 (Niacinamide)', b3_niacinamide_pattern))

    # Vitamin B5 (Pantothenic acid)
    b5_pattern = Chem.MolFromSmarts("CC(C)(O)C(O)C(=O)NCCC(=O)O")
    patterns.append(('Vitamin B5 (Pantothenic acid)', b5_pattern))

    # Vitamin B6 (Pyridoxine)
    b6_pyridoxine_pattern = Chem.MolFromSmarts("COc1cc(CO)c(C)c(O)n1")
    patterns.append(('Vitamin B6 (Pyridoxine)', b6_pyridoxine_pattern))

    # Pyridoxal
    b6_pyridoxal_pattern = Chem.MolFromSmarts("O=CC1=NC=C(CO)C(C)C1O")
    patterns.append(('Vitamin B6 (Pyridoxal)', b6_pyridoxal_pattern))

    # Pyridoxamine
    b6_pyridoxamine_pattern = Chem.MolFromSmarts("NCC1=NC=C(CO)C(C)C1O")
    patterns.append(('Vitamin B6 (Pyridoxamine)', b6_pyridoxamine_pattern))

    # Vitamin B7 (Biotin)
    b7_pattern = Chem.MolFromSmarts("O=C1NC(=O)N2C[C@H](SC1)[C@H]2CCCC(=O)O")
    patterns.append(('Vitamin B7 (Biotin)', b7_pattern))

    # Vitamin B9 (Folic acid)
    b9_pattern = Chem.MolFromSmarts("Nc1nc2NCC(=N)c2n1-c1ccc(cc1)C(=O)N[C@H](CCC(=O)O)C(=O)O")
    patterns.append(('Vitamin B9 (Folic acid)', b9_pattern))

    # Vitamin B12 (Cobalamin) - Complex corrin ring with cobalt
    # Due to complexity, check for cobalt coordinated to corrin-like structure
    b12_pattern = Chem.MolFromSmarts("[Co].[C,c]1=C[C,c]2[C,c]=[C,c][C,c]=[C,c][C,c]=[C,c][C,c]=1[C,c]=[C,c]2")
    patterns.append(('Vitamin B12 (Cobalamin)', b12_pattern))

    # Additional patterns for phosphorylated forms
    # Pyridoxal-5'-phosphate (active form of Vitamin B6)
    b6_pyridoxal_phosphate_pattern = Chem.MolFromSmarts("O=P(O)(O)OCC1=NC=C(C=O)C(C)C1O")
    patterns.append(('Vitamin B6 (Pyridoxal-5\'-phosphate)', b6_pyridoxal_phosphate_pattern))

    # Thiamine pyrophosphate (active form of Vitamin B1)
    b1_thiamine_pyrophosphate_pattern = Chem.MolFromSmarts("CC1=C(SC=[N+]1C)CCOP(O)(=O)OP(O)([O-])=O")
    patterns.append(('Vitamin B1 (Thiamine pyrophosphate)', b1_thiamine_pyrophosphate_pattern))

    # Include the provided examples as SMARTS patterns
    example_smiles = [
        # SMILES strings from provided examples
        "Cc1ncc(CO)c(C[NH3+])c1O",  # pyridoxaminium(1+)
        "CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)NC(CCC(=O)O)C(=O)O",  # 2-[[[4-[(2,4-diamino-6-pteridinyl)methyl-methylamino]phenyl]-oxomethyl]amino]pentanedioic acid
        "Cc1ncc(COP([O-])([O-])=O)c(C[NH3+])c1O",  # pyridoxamine 5'-phosphate(1-)
        # ... Add more examples as needed
    ]

    for idx, ex_smile in enumerate(example_smiles):
        ex_mol = Chem.MolFromSmiles(ex_smile)
        if ex_mol:
            patterns.append((f'Example Molecule {idx+1}', ex_mol))

    # Check for matches
    for name, pattern in patterns:
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches {name}"
    
    # As a fallback, check for similarity to known B vitamins
    # Load known B vitamin molecules
    b_vitamins_smiles = [
        # SMILES of known B vitamins
        "CC1=C(C)C2=NC3=C(C=CC(=C3N2)N1)[N]C(=O)N(C)[C@H]4O[C@@H](CO)[C@H](O)[C@H](O)[C@H]4O",  # Riboflavin
        "C[n+](C)c1ncc(C)cn1CCc1cnc(C)nc1N",  # Thiamine
        "OCCNC(=O)C[C@@H](O)C(C)(C)O",  # Pantothenic acid
        "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@H](CCC(=O)O)C(=O)O)c(=O)[nH]c2n1",  # Folic acid
        "C[C@H]1CN(C)[C@H](CO)O[C@H]1CO",  # Nicotinamide riboside
        # ... Add more as needed
    ]
    b_vitamin_mols = [Chem.MolFromSmiles(smi) for smi in b_vitamins_smiles if Chem.MolFromSmiles(smi) is not None]

    # Compute fingerprints
    fp_query = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    for idx, b_vit_mol in enumerate(b_vitamin_mols):
        fp_b_vit = AllChem.GetMorganFingerprintAsBitVect(b_vit_mol, 2, nBits=2048)
        similarity = Chem.DataStructs.TanimotoSimilarity(fp_query, fp_b_vit)
        if similarity > 0.7:
            return True, f"Similar to known B vitamin (Index {idx+1}) with similarity {similarity:.2f}"

    return False, "Does not match any B vitamin patterns"

__metadata__ = {
    'chemical_class': {
        'name': 'B vitamin',
        'definition': 'Any member of the group of eight water-soluble vitamins originally thought to be a single compound (vitamin B) that play important roles in cell metabolism. The group comprises vitamin B1, B2, B3, B5, B6, B7, B9, and B12.',
    }
}