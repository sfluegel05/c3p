"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for amino group (-[NX3H2])
    amino_pattern = Chem.MolFromSmarts("[NX3H2]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # Check for chirality (except glycine, which is non-chiral)
    chirality_found = False
    glycine_chiral_center_pattern = Chem.MolFromSmarts("C(C(=O)O)N")
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() != Chem.CHI_UNSPECIFIED:
            chirality_found = True
            break
    if not chirality_found and not mol.HasSubstructMatch(glycine_chiral_center_pattern):
        return False, "No chirality found, not glycine"

    # Count alpha carbon connections to ensure there is a proper backbone structure
    valid_alpha_carbon = False
    alpha_carbon_pattern = Chem.MolFromSmarts("N[C@?H]C(=O)O")
    if mol.HasSubstructMatch(alpha_carbon_pattern):
        valid_alpha_carbon = True
    if not valid_alpha_carbon:
        return False, "Invalid alpha-amino acid structure"

    return True, "Contains carboxylic acid and amino group with valid alpha carbon; matches a proteinogenic amino acid"

# Example test case
smiles_example = "N[C@@H](CC(N)=O)C(O)=O"  # L-asparagine
result, reason = is_proteinogenic_amino_acid(smiles_example)
print(f"Is proteinogenic amino acid: {result}, Reason: {reason}")