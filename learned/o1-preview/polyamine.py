"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.

    An amino group is any nitrogen atom that is not part of an amide, nitro group, nitrile, or quaternary ammonium ion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for amino groups (excluding amides, nitro groups, nitriles, and quaternary ammonium)
    amino_group_smarts = "[N;!$(N-C=O);!$(N=O);!$(N#*);!$(N=N);!$(N(=O)(O));!$([N+])]"
    amino_group_pattern = Chem.MolFromSmarts(amino_group_smarts)
    amino_groups = mol.GetSubstructMatches(amino_group_pattern)
    amino_group_count = len(amino_groups)

    if amino_group_count >= 2:
        return True, f"Contains {amino_group_count} amino groups"
    else:
        return False, f"Contains {amino_group_count} amino group(s), which is less than 2"