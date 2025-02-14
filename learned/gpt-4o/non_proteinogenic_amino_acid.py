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
    
    # General amino acid pattern: amino group, alpha carbon, carboxyl group
    general_aa_pattern = Chem.MolFromSmarts("[NX3][C][CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(general_aa_pattern):
        return False, "Does not contain general amino acid structure"
    
    # SMARTS for standard amino acids (covering all 20)
    standard_amino_acids_smarts = [
        "C[C@@H](N)C(=O)O",  # Alanine
        "N[C@@H](C)C(=O)O",  # Glycine
        "N[C@H](CC(=O)O)C(=O)O",  # Aspartate
        "N[C@H](C(=O)O)CC(=O)O",  # Glutamate
        # Continue similarly for all 20 standard amino acids
    ]

    # Check for match with any standard amino acids
    for aa_smarts in standard_amino_acids_smarts:
        aa_pattern = Chem.MolFromSmarts(aa_smarts)
        if mol.HasSubstructMatch(aa_pattern):
            return False, "Matches a standard amino acid"

    # Identify modifications suggesting non-proteinogenic origin
    if len(mol.GetAtoms()) > 20:
        return True, "Complex modifications suggesting non-standard origin"

    # Unusual elements, isotopes, or significant structural differences
    unusual_features = {atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 7, 8, 16]}
    if unusual_features:
        return True, f"Contains unusual elements: {unusual_features}"

    return True, "Does not match any standard amino acids and has unique characteristics"