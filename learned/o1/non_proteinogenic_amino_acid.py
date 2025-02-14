"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI:64049 non-proteinogenic amino acid
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is any amino acid that is not among the standard amino acids
    encoded by the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for amino and carboxylic acid groups
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$([N][C,S]=O)]")  # primary amine not attached to carbonyl/sulfonyl
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-1]")      # carboxylic acid group

    # Check for amino group
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No primary amino group found"

    # Check for carboxylic acid group
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Now, check if the molecule is one of the standard amino acids
    # List of SMILES strings for the 22 standard amino acids (including selenocysteine and pyrrolysine)
    standard_amino_acids_smiles = [
        "N[C@@H](C)C(=O)O",                       # L-Alanine
        "NC(=O)C[C@H](N)C(=O)O",                  # L-Asparagine
        "N[C@@H](CC(=O)O)C(=O)O",                 # L-Aspartic acid
        "N[C@@H](CS)C(=O)O",                      # L-Cysteine
        "NCC(=O)O",                               # Glycine
        "N[C@@H](Cc1c[nH]cn1)C(=O)O",             # L-Histidine
        "N[C@@H](CCCCN)C(=O)O",                   # L-Lysine
        "N[C@@H](CCC(N)=O)C(=O)O",                # L-Glutamine
        "N[C@@H](CCC(=O)O)C(=O)O",                # L-Glutamic acid
        "N[C@@H](Cc1ccccc1)C(=O)O",               # L-Phenylalanine
        "N[C@@H](CO)C(=O)O",                      # L-Serine
        "N[C@@H](C(C)O)C(=O)O",                   # L-Threonine
        "N[C@@H](Cc1c[cH][cH]c1)C(=O)O",          # L-Tryptophan
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O",            # L-Tyrosine
        "N[C@@H](CC(C)C)C(=O)O",                  # L-Leucine
        "N[C@@H](C(C)C)C(=O)O",                   # L-Valine
        "N[C@@H](CC(C)CC)C(=O)O",                 # L-Isoleucine
        "N[C@@H](CSCC[SeH])C(=O)O",               # L-Selenocysteine
        "N[C@@H](CC1CN1)C(=O)O",                  # L-Proline
        "N[C@@H](CCCNC(N)=N)C(=O)O",              # L-Arginine
        "N[C@@H](C[S](=O)C)C(=O)O",               # L-Methionine
        "N[C@@H](CC1=CN=CN1)C(=O)O",              # L-Pyrrolysine
    ]

    # Check if the molecule matches any of the standard amino acids
    # Use molecular fingerprints for more robust matching
    mol_fp = AllChem.GetMorganFingerprint(mol, radius=2)
    for std_smiles in standard_amino_acids_smiles:
        std_mol = Chem.MolFromSmiles(std_smiles)
        if std_mol is None:
            continue
        std_fp = AllChem.GetMorganFingerprint(std_mol, radius=2)
        similarity = DataStructs.TanimotoSimilarity(mol_fp, std_fp)
        if similarity == 1.0:
            return False, "Molecule is a standard (proteinogenic) amino acid"

    return True, "Molecule is a non-proteinogenic amino acid (contains amino and carboxyl groups but is not one of the standard amino acids)"