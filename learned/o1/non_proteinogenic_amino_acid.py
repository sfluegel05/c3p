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

    # Define SMARTS patterns for alpha-amino acid backbone
    # Matches an alpha carbon with amino group and carboxyl group
    alpha_amino_acid_pattern = Chem.MolFromSmarts("""
        [C;!$(C=O);!$(C#N);!$(C=[N,O,S])](N)(C(=O)[O,H,-1])
    """)

    # Check for alpha-amino acid backbone
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    if not matches:
        return False, "No alpha-amino acid backbone found"

    # Ensure there's only one alpha-amino acid backbone to exclude peptides
    if len(matches) > 1:
        return False, f"Multiple amino acid backbones found ({len(matches)}), possible peptide"

    # Now, check if the molecule is one of the standard amino acids
    # List of standard amino acids as RDKit molecules
    standard_amino_acids_smiles = [
        "N[C@@H](C)C(=O)O",                       # L-Alanine
        "N[C@@H](CC(=O)O)C(=O)O",                 # L-Aspartic acid
        "N[C@@H](C(=O)N)C(=O)O",                  # L-Asparagine
        "N[C@@H](CS)C(=O)O",                      # L-Cysteine
        "NCC(=O)O",                               # Glycine
        "N[C@@H](CCC(=O)O)C(=O)O",                # L-Glutamic acid
        "N[C@@H](CC(=O)N)C(=O)O",                 # L-Glutamine
        "N[C@@H](Cc1c[nH]cn1)C(=O)O",             # L-Histidine
        "N[C@@H](I)C(=O)O",                       # L-Isoleucine
        "N[C@@H](CCCN)C(=O)O",                    # L-Lysine
        "N1CCCC1C(=O)O",                          # L-Proline
        "N[C@@H](C(C)C)C(=O)O",                   # L-Valine
        "N[C@@H](CC1=CNC=N1)C(=O)O",              # L-Tryptophan
        "N[C@@H](Cc1ccccc1)C(=O)O",               # L-Phenylalanine
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O",            # L-Tyrosine
        "N[C@@H](CSCC[SeH])C(=O)O",               # L-Selenocysteine
        "N[C@@H](C(=O)O)C(=O)O",                  # L-Aspartic acid (beta)
        # ... (add other standard amino acids)
    ]

    standard_amino_acids = [Chem.MolFromSmiles(smi) for smi in standard_amino_acids_smiles]

    # Check for exact match with standard amino acids (including stereochemistry)
    for std_mol in standard_amino_acids:
        if std_mol is None:
            continue
        if mol.HasSubstructMatch(std_mol):
            return False, "Molecule is a standard (proteinogenic) amino acid"

    # Check for peptide bonds to exclude peptides
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    # Subtract 1 because the amino acid itself will have one amide bond if it's an amide
    if len(peptide_bonds) > 1:
        return False, "Peptide bonds detected, molecule may be a peptide"

    return True, "Molecule is a non-proteinogenic amino acid (contains alpha-amino acid backbone but is not a standard amino acid)"