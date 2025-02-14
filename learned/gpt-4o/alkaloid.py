"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determine if a molecule is likely an alkaloid based on its SMILES string.
    This function attempts to identify complex heterocyclic structures with nitrogen,
    focusing on the diversity of potential alkaloid structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define broader patterns for known alkaloid-like nitrogen-containing ring systems
    alkaloid_substructures = [
        Chem.MolFromSmarts('c1ncccc1'),  # Pyridine-like
        Chem.MolFromSmarts('n1ccccc1'),  # Pyrrole-like
        Chem.MolFromSmarts('Cn1cccc1'),  # Tertiary nitrogen in aromatic
        Chem.MolFromSmarts('C1CNCCN1'),  # Piperidine/pyrrolidine-like
        Chem.MolFromSmarts('c1[nH]c2ccccc2c1'),  # Indole-like structures
        Chem.MolFromSmarts('C1CCN(C1)C2=CC=CC=C2'),  # Quinolizidine-like
        Chem.MolFromSmarts('C1CN2CCC1CC2'),  # Tropane-like
        Chem.MolFromSmarts('c1c2[nH]c3ccccc3c2[nH]c1'), # Isoquinoline-like
        Chem.MolFromSmarts('c1[nH]c2ncccc2c1')  # Purine-like
    ]

    # Check for nitrogen in recognized alkaloid substructures
    for idx, substructure in enumerate(alkaloid_substructures):
        if mol.HasSubstructMatch(substructure):
            return True, f"Contains recognized alkaloid-like nitrogen ring system #{idx + 1}"

    # Exclude simple amines or peptide-like structures
    peptide_pattern = Chem.MolFromSmarts('C(=O)N-C')
    amine_pattern = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')  # Primary and secondary amines
    if mol.HasSubstructMatch(peptide_pattern) or mol.HasSubstructMatch(amine_pattern):
        return False, "Contains simple amine or peptide-like structure, likely not an alkaloid"

    return False, "No recognizable alkaloid structure found"

# Example usage:
# smiles = "CCCCCCCCCCC1=CC=C(N1)\\C=C1/N=C(C=C1OC)C1=CC=CN1"
# result, reason = is_alkaloid(smiles)
# print(result, reason)