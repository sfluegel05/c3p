"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine consists of a sphingosine backbone with an N-linked acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for sphingosine backbone: part of long chain with alcohols and double bond structure
    sphingosine_pattern = Chem.MolFromSmarts("[NH][C@@H]([CH2]O)[C@H](O)C=C") 
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Pattern for N-acyl group: N-attached acyl group (NH-C(=O)-C)
    n_acyl_pattern = Chem.MolFromSmarts("N[C;X4](=O)[C;X4]")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group found"
    
    return True, "Contains sphingosine backbone with N-linked acyl group"

# Example usage
smiles_example = "CCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
result, reason = is_N_acylsphingosine(smiles_example)
print(result, reason)