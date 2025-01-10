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

    # Define an expanded pattern for D-ribose, allowing modifications like methoxy
    ribose_pattern = Chem.MolFromSmarts("O[c@H]1[c@H](O)[c@@H](C)[c@@H](O)[C@H](O)1")  # More flexible ribose pattern
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar component found"

    # Define broader patterns for nucleobases in ribonucleosides
    extended_nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2n(cnc12)"),  # Purine base pattern
        Chem.MolFromSmarts("c1cnc[nH]c1"),     # Pyrimidine base pattern
        Chem.MolFromSmarts("c1c[nH]c(=O)n(c1)"),  # Modified pyrimidine pattern
        Chem.MolFromSmarts("n1cnc2c(ncnc2)n1"), # Include potential alternate purine patterns
        Chem.MolFromSmarts("n1cncc1N"),        # Generalized pattern to include possible modifications
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