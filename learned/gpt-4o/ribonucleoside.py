"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is any nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a pattern for D-ribose, allowing basic chemical modifications, like hydroxy or methyl groups
    ribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@H](CO)1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar component found"

    # Define broader patterns for nucleobases
    extended_nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2n(cnc12)"),  # Purine
        Chem.MolFromSmarts("c1cnc[nH]c1"),     # Pyrimidine
        Chem.MolFromSmarts("c1c[nH]c(=O)n(c1)"), # Modified pyrimidine
    ]

    # Check for nucleobase attachment
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in extended_nucleobase_patterns)

    if has_nucleobase:
        return True, "Contains D-ribose sugar with an attached nucleobase"
    else:
        return False, "No nucleobase found attached to the ribose"

# Example usage
# smiles = "CN(C)c1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](N)[C@H]1O"  # Example ribonucleoside SMILES
# result, reason = is_ribonucleoside(smiles)
# print(result, reason)