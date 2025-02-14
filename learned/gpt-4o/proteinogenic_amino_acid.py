"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    Checks for amino acids with a glycine structure or standard L-alpha-amino acids with
    standalone carboxyl and amino groups.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES to create RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycine special case: it's non-chiral
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    if mol.HasSubstructMatch(glycine_pattern):
        return True, "Glycine detected, which is an achiral proteinogenic amino acid"

    # Look for L-alpha-amino acid pattern including carboxyl and amino groups on alpha carbon
    # and ensuring the alpha carbon is chiral or matches common configurations for L-amino acids.
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[N;!H0][C@@H](C(=O)O)[C,!N]")
    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        # Check for modifications like selenocysteine or pyrrolysine
        modifications = ["Se", "Pyrrolysine", "formyl"]
        for mod in modifications:
            if mod in smiles:
                return True, f"Detected proteinogenic amino acid modification: {mod}"
        return True, "L-alpha-amino acid detected, a standard proteinogenic amino acid"

    # Exclude peptides and modified amino acids by checking entity size and typical bonds
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N[C]")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Peptide bond detected, molecule may be part of a peptide"

    return False, "Not a standalone proteinogenic amino acid"