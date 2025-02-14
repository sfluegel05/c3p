"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least two distinct biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a bioconjugate, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define biological components using SMARTS patterns (simplified for the example)
    # Amino acid backbone pattern (e.g., C-C(C(=O)O)N)
    amino_acid_pattern = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")
    
    # Fatty acid chain pattern (e.g., long aliphatic chain with carboxyl group)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)[OX1][C;R0]!@[C;R0]!@[C;R0]")

    # Nucleotide pattern (simplified sugar base)
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")

    # Identify presence of these components
    has_amino_acid = mol.HasSubstructMatch(amino_acid_pattern)
    has_fatty_acid = mol.HasSubstructMatch(fatty_acid_pattern)
    has_nucleotide = mol.HasSubstructMatch(nucleotide_pattern)

    # Check for at least two distinct biological components
    components = [has_amino_acid, has_fatty_acid, has_nucleotide]
    if sum(components) >= 2:
        return True, "Contains at least two different biological moieties covalently linked"
    else:
        return False, "Does not contain at least two distinct biological components"

# Example use (replace SMILES with an actual string to test):
# result, reason = is_bioconjugate("your_smiles_here")