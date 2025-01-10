"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.

    An alpha-amino acid ester has: 
    - An alpha-carbon connected to an amino group and an esterified carboxylic acid group.

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

    # Define SMARTS patterns for alpha-amino acid ester components
    # Allow flexibility in recognizing amine types (primary, secondary, tertiary)
    amino_acid_ester_pattern = Chem.MolFromSmarts('[CX4H1]([NX3])-[CX3](=O)-[OX2][CX4]')
    if not mol.HasSubstructMatch(amino_acid_ester_pattern):
        return False, "No alpha-amino acid ester pattern found"

    # Check for presence of an ester linkage [-C(=O)-O-]
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"

    # Check for presence of a nitrogen atom attached to the alpha-carbon
    # Ensure nitrogen is part of an amine (but can be primary, secondary, or tertiary)
    amino_group_pattern = Chem.MolFromSmarts('[NX3]')
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No amino group found"

    return True, "Molecule contains an alpha-amino acid ester structure"

# Example usage
example_smiles = "COC(=O)CN"
result, reason = is_alpha_amino_acid_ester(example_smiles)
print(f"Classification result: {result}, Reason: {reason}")