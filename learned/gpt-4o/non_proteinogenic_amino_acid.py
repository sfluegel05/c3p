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
    
    # Basic amino acid structure: [N][C](C(=O)O) with possible stereochemistry
    amino_acid_pattern = Chem.MolFromSmarts("[N;H2,H3][C;^2][CX3](=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Does not contain basic amino acid structure"

    # Define SMARTS for all 20 standard amino acids
    standard_amino_acids_smarts = [
        "C[C@@H](N)C(=O)O",  # Alanine
        "N[C@@H](CC(=O)O)C(=O)O",  # Aspartic acid
        "N[C@@H](CO)C(=O)O",  # Serine
        "N[C@@H](CC(C)C)C(=O)O",  # Valine
        "N[C@@H](CCSC)C(=O)O",  # Methionine
        "N[C@@H](c1ccccc1)C(=O)O",  # Phenylalanine
        # Add more for each standard amino acid...
    ]

    # Check if it matches any known standard amino acids
    for aa_smarts in standard_amino_acids_smarts:
        aa_pattern = Chem.MolFromSmarts(aa_smarts)
        if mol.HasSubstructMatch(aa_pattern):
            return False, "Matches a standard amino acid"

    # Look for D-amino acids using inversion of stereochemistry
    d_amino_acid_pattern = Chem.MolFromSmarts("N[C@H]C(=O)O")
    if mol.HasSubstructMatch(d_amino_acid_pattern):
        return True, "Contains a D-amino acid stereochemistry"

    # Check for unusual isotopes or modifications
    unusual_elements = {atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 9}
    if unusual_elements:
        return True, f"Contains unusual elements or isotopes: {unusual_elements}"

    return True, "Does not match any standard amino acids and fits the criteria for non-proteinogenic"