"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester typically has the O=C(O)-O- functional group, where
    the central carbon is doubly bonded with an oxygen and singly bonded with
    an alkoxy or aryloxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to get molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # A more precise set of SMARTS patterns that capture carbonate ester core
    patterns = [
        Chem.MolFromSmarts("[$([O][C](=[O])[O][C])]"), # Linear carbonate ester O=C(O)O-C
        Chem.MolFromSmarts("[$([O][C](=[O])[O][CX3])]"), # Cyclic carbonate, considering sp3 carbons
        Chem.MolFromSmarts("C(=O)O[C;H1]"), # Ensure one oxygen has at least one hydrogen or a direct carbon bond for diversity
        Chem.MolFromSmarts("O=C1OC1"), # Base small cyclic carbonate
    ]

    # Check if the molecule matches any pattern
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbonate ester functional group"
    
    # Evaluate failure reasons by analyzing functional groups on unusual positives
    c_d = Chem.MolFromSmarts("[$([C]=[O])]")
    o_d = Chem.MolFromSmarts("[#8]") # Standalone [O], if there is ambiguity
    if mol.HasSubstructMatch(c_d) and mol.HasSubstructMatch(o_d):
        return False, "Contains C=O and O groups but not in Typical Ester Group"

    return False, "No carbonate ester functional group found"