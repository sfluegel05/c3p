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

    # Define SMARTS pattern for alpha-amino acid backbone
    # Matches an alpha carbon with amino group and carboxyl group
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2][C@@H]([*])C(=O)[O-]")
    if alpha_amino_acid_pattern is None:
        return False, "Invalid SMARTS pattern for alpha-amino acid backbone"

    # Add hydrogens (necessary for accurate matching)
    mol = Chem.AddHs(mol)

    # Check for alpha-amino acid backbone
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    if not matches:
        return False, "No alpha-amino acid backbone found"

    # Ensure there's only one alpha-amino acid backbone to exclude peptides
    if len(matches) > 1:
        return False, f"Multiple amino acid backbones found ({len(matches)}), possible peptide"

    # Generate canonical SMILES for the molecule
    mol_canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    # List of standard amino acids canonical SMILES (with stereochemistry)
    standard_amino_acids_smiles = [
        "N[C@@H](C)C(=O)O",                          # L-Alanine
        "N[C@@H](CC(=O)O)C(=O)O",                    # L-Aspartic acid
        "N[C@@H](C(=O)N)C(=O)O",                     # L-Asparagine
        "N[C@@H](CS)C(=O)O",                         # L-Cysteine
        "NCC(=O)O",                                  # Glycine
        "N[C@@H](CCC(=O)O)C(=O)O",                   # L-Glutamic acid
        "N[C@@H](CC(=O)N)C(=O)O",                    # L-Glutamine
        "N[C@@H](Cc1c[nH]cn1)C(=O)O",                # L-Histidine
        "N[C@@H](I)C(=O)O",                          # L-Isoleucine
        "N1CCCC1C(=O)O",                             # L-Proline
        "N[C@@H](CC(C)C)C(=O)O",                     # L-Leucine
        "N[C@@H](CO)C(=O)O",                         # L-Serine
        "N[C@@H](Cc1ccccc1)C(=O)O",                  # L-Phenylalanine
        "N[C@@H](Cc1ccc(O)cc1)C(=O)O",               # L-Tyrosine
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",          # L-Tryptophan
        "N[C@@H](CCCN)C(=O)O",                       # L-Lysine
        "N[C@@H](CCS)C(=O)O",                        # L-Methionine
        "N[C@@H](C(C)O)C(=O)O",                      # L-Threonine
        "N[C@@H](C(C)C)C(=O)O",                      # L-Valine
        "N[C@@H](C(=O)O)C(=O)O",                     # L-Serine (alternative representation)
        # Add any missing standard amino acids
    ]

    # Generate canonical SMILES for standard amino acids
    standard_amino_acids_canon_smiles = set()
    for smi in standard_amino_acids_smiles:
        std_mol = Chem.MolFromSmiles(smi)
        if std_mol:
            std_mol_canon = Chem.MolToSmiles(std_mol, isomericSmiles=True)
            standard_amino_acids_canon_smiles.add(std_mol_canon)

    # Check if the molecule matches any standard amino acid
    if mol_canonical_smiles in standard_amino_acids_canon_smiles:
        return False, "Molecule is a standard (proteinogenic) amino acid"

    # Check for peptide bonds to exclude peptides
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N[C;!$([C;!H0])]")
    if peptide_bond_pattern is None:
        return False, "Invalid SMARTS pattern for peptide bond"

    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if peptide_bonds:
        return False, "Peptide bonds detected, molecule may be a peptide"

    return True, "Molecule is a non-proteinogenic amino acid (contains alpha-amino acid backbone but is not a standard amino acid)"