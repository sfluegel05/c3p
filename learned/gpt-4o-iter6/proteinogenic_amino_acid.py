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

    # Standard amino acid core patterns
    amino_acid_pattern = Chem.MolFromSmarts("[N+0;!H2][C@&H1;X4]([C;X3](=[O;X1])[O;X2])[!#1]")
    glycine_pattern = Chem.MolFromSmarts("[N+0]CC(=O)[O;X2]")

    # Check if the mol has the amino acid backbone (allows isotopic labels/deuteration)
    if not mol.HasSubstructMatch(amino_acid_pattern):
        if mol.HasSubstructMatch(glycine_pattern):
            return True, "SMILES matches glycine structure, an exception to chirality"
        return False, "No standard amino acid backbone with valid chirality found"

    # Further chirality and derivative check
    chiral_atoms = [atom for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.CHI_UNSPECIFIED]
    if not chiral_atoms:
        return False, "Lacking required chirality or has unresolved chirality"

    # Ensure single amino acid and check for lone peptide bonds or extended structures
    if mol.GetNumBonds() > 10:  # Arbitrary threshold for complexity beyond single amino acid
        return False, "Molecule too complex; likely a peptide or larger compound"

    return True, "Valid backbone and chirality detected; matches a proteinogenic amino acid"

# Example test case
smiles_example = "N[C@@H](CC(N)=O)C(O)=O"  # L-asparagine
result, reason = is_proteinogenic_amino_acid(smiles_example)
print(f"Is proteinogenic amino acid: {result}, Reason: {reason}")