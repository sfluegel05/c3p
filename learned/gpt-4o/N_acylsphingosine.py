"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: N-acylsphingosine
"""
from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is composed of sphingosine having a fatty acyl group attached to the nitrogen.

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

    # Look for sphingosine backbone pattern:
    # Simplified pattern considering mostly the backbone
    sphingosine_pattern = Chem.MolFromSmarts("O[C@@H](CO)C[C@@H](C=C)CCCC(=N)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Check for presence of double bond in the backbone
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "Missing double bond in sphingosine-like structure"
    
    # Look for amide bond pattern where nitrogen is connected to carbonyl carbon of a fatty acyl chain
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found connecting nitrogen to fatty acyl group"

    return True, "Contains sphingosine backbone with an amide bond connecting to a fatty acyl group"