"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid typically contains an indole skeleton, modified by various alkaloidal features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More relaxed pattern to search for variations of the indole skeleton
    generic_indole_pattern = Chem.MolFromSmarts('c1c[nH]c2cccc2c1')  # Allow more variations in the indole ring
    extended_indole_pattern = Chem.MolFromSmarts('c1cnc2cccc2c1')

    # Check if either general indole or expanded pattern is present
    if not (mol.HasSubstructMatch(generic_indole_pattern) or mol.HasSubstructMatch(extended_indole_pattern)):
        return False, "No recognizable indole-like skeleton found"

    # Check for at least one additional nitrogen, beyond what indole naturally contains
    additional_nitrogen_pattern = Chem.MolFromSmarts('[#7]')  # Any nitrogen atom
    additional_n_count = len(mol.GetSubstructMatches(additional_nitrogen_pattern))

    # Typically, for indole alkaloids, more than one nitrogen is common
    if additional_n_count < 2:
        return False, "Insufficient additional nitrogen atoms, found only indole-based nitrogen"

    return True, "Molecule possesses an indole-like structure with characteristics of an alkaloid"

# Example usage for testing based on provided SMILES strings
example_smiles = "C1=C2CC[NH2+]3C4=CC=CC=C4NC2=CC3=NC1C=C"
result, reason = is_indole_alkaloid(example_smiles)
print(result, reason)