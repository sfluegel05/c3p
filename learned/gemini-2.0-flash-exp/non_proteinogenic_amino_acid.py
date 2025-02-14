"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is an amino acid not naturally encoded in the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Basic amino acid structure check: N-C-C(=O)-O
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3;H2,H1][CX4;H1]([H])[CX3](=[OX1])[OX2;H1,H0]")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "Missing basic amino acid structure"

    # 2. Check if it is a standard proteinogenic amino acid
    # Use SMARTS to match each of the standard amino acid sidechains
    # Explicitly account for L and D forms (where appropriate).
    proteinogenic_amino_acids_smarts = [
        "[NX3;H2,H1][C@H]([CH3])[CX3](=[OX1])[OX2;H1,H0]",  # L-Alanine
        "[NX3;H2,H1][C@@H](CC(=O)N)[CX3](=[OX1])[OX2;H1,H0]",  # L-Asparagine
        "[NX3;H2,H1][C@H](CC(=O)O)[CX3](=[OX1])[OX2;H1,H0]",  # L-Aspartic acid
        "[NX3;H2,H1][C@@H](CC(=O)O)[CX3](=[OX1])[OX2;H1,H0]",  # D-Aspartic acid
        "[NX3;H2,H1][C@H](Cc1c[nH]cn1)[CX3](=[OX1])[OX2;H1,H0]",  # L-Histidine
        "[NX3;H2,H1][C@@H](Cc1ccccc1)[CX3](=[OX1])[OX2;H1,H0]",  # L-Phenylalanine
        "[NX3;H2,H1][C@H](CS)[CX3](=[OX1])[OX2;H1,H0]",  # L-Cysteine
        "[NX3;H2,H1][C@@H](CC(=O)N)[CX3](=[OX1])[OX2;H1,H0]",  # L-Glutamine
        "[NX3;H2,H1][C@H](CCC(=O)O)[CX3](=[OX1])[OX2;H1,H0]",  # L-Glutamic acid
        "[NX3;H2,H1][C@@H](CCC(=O)O)[CX3](=[OX1])[OX2;H1,H0]",  # D-Glutamic acid
        "[NX3;H2,H1][C@H](C(C)C)[CX3](=[OX1])[OX2;H1,H0]",  # L-Valine
        "[NX3;H2,H1][C@H](CCCNC(=N)N)[CX3](=[OX1])[OX2;H1,H0]", # L-Arginine
        "[NX3;H2,H1][C@H](CO)[CX3](=[OX1])[OX2;H1,H0]",  # L-Serine
        "[NX3;H2,H1][C@H]([C@@H](O)C)[CX3](=[OX1])[OX2;H1,H0]",  # L-Threonine
        "[NX3;H2,H1][C@H](CCSC)[CX3](=[OX1])[OX2;H1,H0]",  # L-Methionine
        "[NX3;H2,H1][C@H](C(C)CC)[CX3](=[OX1])[OX2;H1,H0]",  # L-Leucine
        "[NX3;H2,H1][C@H](C(C)(C)C)[CX3](=[OX1])[OX2;H1,H0]",  # L-Isoleucine
        "[NX3;H2,H1][C@H](Cc1c[nH]c2ccccc12)[CX3](=[OX1])[OX2;H1,H0]", # L-Tryptophan
        "[NX3;H2,H1][C@H](Cc1ccc(O)cc1)[CX3](=[OX1])[OX2;H1,H0]",  # L-Tyrosine
        "N1[C@@H]2[C@H](C[C@@H]1C(=O)O)N[CH2]2",  # L-Proline
    ]

    for smarts in proteinogenic_amino_acids_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            return False, "It is a standard proteinogenic amino acid."

    return True, "It is a non-proteinogenic amino acid."