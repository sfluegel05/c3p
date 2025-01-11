"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.

    An alpha-amino acid ester typically has:
    - An alpha-carbon connected to an amino group.
    - A carboxyl group (-C(=O)O-) esterified.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern: capture alpha-amino acid ester structures
    # Pattern for alpha-amino acid ester [N]C(C(=O)O)[CX4,CX3]
    alpha_amino_ester_patterns = [
        Chem.MolFromSmarts('[NX3][CX4;H0,H1](C(=O)[OX2][C])'),  # Linear form
        Chem.MolFromSmarts('C(=O)[OX2][CX4][NX3][CX4]'),  # Potentially cyclic including aromatics
        Chem.MolFromSmarts('[C,N]C([NX3])[CX4,CX3](=O)[O][C]'),  # More generic pattern
    ]

    # Check for each pattern
    for pattern in alpha_amino_ester_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains an alpha-amino acid ester structure"
    
    return False, "No alpha-amino acid ester pattern found"

# Example usage
example_smiles = "COC(=O)CN"
result, reason = is_alpha_amino_acid_ester(example_smiles)
print(f"Classification result: {result}, Reason: {reason}")