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

    # Look for standard amino acid backbone (N-C-C(=O)-O)
    backbone_pattern = Chem.MolFromSmarts("N[C@?H](C(=O)O)")

    # Check if the mol has the standard amino acid backbone
    if not mol.HasSubstructMatch(backbone_pattern):
        # Glycine exception, which is achiral
        glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
        if mol.HasSubstructMatch(glycine_pattern):
            return True, "Matches glycine structure, an exception to chirality"
        return False, "No standard amino acid backbone structure found"

    # Check for chirality (except glycine, which is non-chiral)
    chirality_found = False
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() != Chem.CHI_UNSPECIFIED:
            chirality_found = True
            break

    if not chirality_found:
        return False, "No chirality found except for glycine"

    return True, "Contains valid backbone and chirality; matches a proteinogenic amino acid"

# Example test case
smiles_example = "N[C@@H](CC(N)=O)C(O)=O"  # L-asparagine
result, reason = is_proteinogenic_amino_acid(smiles_example)
print(f"Is proteinogenic amino acid: {result}, Reason: {reason}")