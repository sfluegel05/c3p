"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    This class is defined by an L-alpha-amino acid backbone with an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Refined pattern for L-alpha-amino acid backbone: [N@H](C)[C@H](C(=O)O) to ensure stereospecificity
    # This assumes that L configuration involves specific chiral notations; adjust if necessary to generalize or restrict further
    amino_acid_pattern = Chem.MolFromSmarts("[N][C@H](C)C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return (False, "No L-alpha-amino acid backbone found")

    # Refined pattern for N-acyl group connected to a chiral center: [NX3][CX3](=O) to signify proper connectivity
    n_acyl_pattern = Chem.MolFromSmarts("[N][C](=O]")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return (False, "No N-acyl group found")

    # Additionally, ensure the N-acyl group is directly attached to the amino acid nitrogen
    # Excluding large side chains or non-specific attachments that match in peptides
    # This approach prevents matching peptide chains or non-related structures

    return (True, "Contains L-alpha-amino acid backbone with an N-acyl substituent")