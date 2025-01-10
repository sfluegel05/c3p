"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine has a phytosphingosine backbone with a fatty acyl group attached to nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General pattern for a phytosphingosine backbone
    # Long chain with a primary amino group (N) and at least two hydroxyl (OH) groups
    phytosphingosine_pattern = Chem.MolFromSmarts("C[C@H](O)C[C@H](O)CN")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Look for fatty acyl group: amide linkage to the nitrogen
    # R-C(=O)-N, where R is a long chain (CX4 denotes sp3 carbon)
    fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)N[C]")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No fatty acyl group attached to nitrogen found"

    return True, "Contains phytosphingosine backbone with fatty acyl group attached to nitrogen"