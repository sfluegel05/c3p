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

    # Specific pattern for amino acids (capturing standard structure including isotopically labeled versions)
    standard_amino_acid_pattern = Chem.MolFromSmarts("[N;!H2][C@@,C@&!H0][C;!H0](=O)[O;!H0]")
    strict_glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")

    # Check for basic structure match, considering cases for glycine
    if mol.HasSubstructMatch(standard_amino_acid_pattern):
        # Further evaluate complexity of the molecule
        if any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms() if mol.GetBondsBetweenAtoms(atom.GetIdx(), atom.GetIdx())):
            return False, "Contains potentially extended carbon backbone beyond amino acid core"
        return True, "Valid amino acid backbone and chirality detected"
    elif mol.HasSubstructMatch(strict_glycine_pattern):
        return True, "SMILES matches glycine structure, the exception for non-chirality"
    
    return False, "Does not conform to established proteinogenic amino acid structure"

# Example test case
smiles_example = "N[C@@H](CC(N)=O)C(O)=O"  # L-asparagine
result, reason = is_proteinogenic_amino_acid(smiles_example)
print(f"Is proteinogenic amino acid: {result}, Reason: {reason}")