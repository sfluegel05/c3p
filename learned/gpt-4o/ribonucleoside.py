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
    
    # Identify the ribose sugar pattern
    ribose_pattern = Chem.MolFromSmarts("C1[C@@H](O)[C@@H](O)[C@H](CO)O1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar component found"
    
    # Identify a common nucleobase attached (for simplicity, we use a generic pattern)
    nucleobase_pattern = Chem.MolFromSmarts("n1cnc2c1[nH]c(=O)[nH]c2=O")
    purine_base_pattern = Chem.MolFromSmarts("c1ncnc2n1cnc2")
    pyrimidine_base_pattern = Chem.MolFromSmarts("c1[nH]c(=O)nc(=O)[nH]c1")

    has_nucleobase = mol.HasSubstructMatch(nucleobase_pattern) or \
                     mol.HasSubstructMatch(purine_base_pattern) or \
                     mol.HasSubstructMatch(pyrimidine_base_pattern)

    if not has_nucleobase:
        return False, "No nucleobase found attached to the ribose"

    return True, "Contains D-ribose sugar with an attached nucleobase"

# Example usage
# smiles = "CN(C)c1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](N)[C@H]1O"  # Example ribonucleoside SMILES
# result, reason = is_ribonucleoside(smiles)
# print(result, reason)