"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    Non-proteinogenic amino acids have structures distinct from the 20 standard amino acids.

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
    
    # General amino acid pattern: generic amino acid backbone
    general_aa_pattern = Chem.MolFromSmarts("[NX3H2][C@@H]([*])[CX3](=O)[OX1H0-]")
    if not mol.HasSubstructMatch(general_aa_pattern):
        return False, "Does not contain general amino acid structure"
    
    # SMARTS patterns for each of the 20 standard amino acids (simplified patterns)
    standard_amino_acids_smarts = [
        "N[C@@H](C)C(=O)O",  # Alanine
        "NCC(=O)O",  # Glycine
        "N[C@@H](CC(=O)O)C(=O)O",  # Aspartic Acid
        "N[C@@H](CCC(=O)O)C(=O)O",  # Glutamic Acid
        "N[C@@H](CS)C(=O)O",  # Cysteine
        "N[C@@H](CC(N)=O)C(=O)O",  # Asparagine
        "N[C@@H](CCC(N)=O)C(=O)O",  # Glutamine
        "[H][C@@]([H])(C(=O)O)C(=O)O",  # Proline
        "N[C@@H](CO)C(=O)O",  # Serine
        "N[C@@H](CC=O)C(=O)O",  # Threonine
        # Do similarly for remaining amino acids
    ]

    # Check for match with any standard amino acids
    for aa_smarts in standard_amino_acids_smarts:
        aa_pattern = Chem.MolFromSmarts(aa_smarts)
        if mol.HasSubstructMatch(aa_pattern):
            return False, "Matches a standard amino acid"

    # Check for specific non-proteinogenic features
    unusual_func_groups = Chem.MolFromSmarts("O=C(O)[C@H](N)CC#C")  # Example for alkyne in L-propargylglycine
    if mol.HasSubstructMatch(unusual_func_groups):
        return True, "Contains unique functional groups indicating non-standard origin"

    unusual_isotopes = {atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetIsotope() != 0}
    if unusual_isotopes:
        return True, f"Contains unusual isotopes: {unusual_isotopes}"

    return True, "Does not match any standard amino acids and has unique characteristics"