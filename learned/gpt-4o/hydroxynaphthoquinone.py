"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is defined as any naphthoquinone in which the naphthoquinone moiety
    is substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more generalized SMARTS pattern for naphthoquinone core
    # Naphthoquinone structure with flexibility for additional aromatic substitutions
    naphthoquinone_pattern = Chem.MolFromSmarts("C1=CC=C2C(=C1)C=CC(=O)C2=O")
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No naphthoquinone core structure found"
    
    # Define SMARTS pattern for hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy group found"

    return True, "Contains naphthoquinone core with at least one hydroxy group attached"

# Example Usage:
# smiles_example = "O=C1C=CC(=O)c2cc(O)ccc12"  # Example of a hydroxynaphthoquinone
# result, reason = is_hydroxynaphthoquinone(smiles_example)
# print(result, reason)