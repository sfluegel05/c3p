"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool, str: True and reason if molecule is a proteinogenic amino acid, otherwise False and reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define patterns for amino acid backbones
    chiral_amino_acid_pattern = Chem.MolFromSmarts("[N;!H2][C@;R][C;!H0](=O)[O;!H0]")
    achiral_glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    
    # Check for standard chiral amino acids
    if mol.HasSubstructMatch(chiral_amino_acid_pattern):
        return True, "Standard chiral amino acid detected"
    
    # Check for glycine, which is achiral
    if mol.HasSubstructMatch(achiral_glycine_pattern):
        return True, "Glycine detected"

    return False, "No matching pattern for proteinogenic amino acids"

# Example test case
smiles_example = "N[C@@H](CC(N)=O)C(O)=O"  # L-asparagine
result, reason = is_proteinogenic_amino_acid(smiles_example)
print(f"Is proteinogenic amino acid: {result}, Reason: {reason}")