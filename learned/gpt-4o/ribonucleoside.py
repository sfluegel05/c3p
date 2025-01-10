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
    
    # Identify a more flexible ribose sugar pattern allowing minor modifications
    ribose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@H](CO)[C@@H](O)O1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar component found"

    # Patterns for extended nucleobase matching, accounting for common modifications
    extended_nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2n1cnc2"),  # Purines 
        Chem.MolFromSmarts("c1[nH]c(=O)nc(=O)[nH]c1"),  # Pyrimidines 
        Chem.MolFromSmarts("n1cnc2c1[nH]c(=O)[nH]c2=O"),  # Alternative purine/pyrimidine structures
        Chem.MolFromSmarts("c1c[nH]n[cH]1")  # Allow for modifications
    ]

    # Check for any of the nucleobase patterns
    for pattern in extended_nucleobase_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains D-ribose sugar with an attached nucleobase"

    return False, "No nucleobase found attached to the ribose"

# Example usage
# smiles = "CN(C)c1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](N)[C@H]1O"  # Example ribonucleoside SMILES
# result, reason = is_ribonucleoside(smiles)
# print(result, reason)