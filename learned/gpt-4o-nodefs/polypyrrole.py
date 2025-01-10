"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A true polypyrrole is defined as a conjugated polymeric chain of pyrrole units linked via carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for polymeric pyrrole chain (several pyrroles linked in a linear or branched fashion)
    polypyrrole_pattern = Chem.MolFromSmarts('n1cccc1-!@[C]-n2cccc2')
    
    # Cyclic structure exclusion to avoid porphyrin-like structures
    cyclic_pyrrole_exclusion = Chem.MolFromSmarts('[$(c1cc[nH]c[c2ncc3nc[nH]c3[nH]2]c1)]')
    
    # Check if any pattern relevant to polypyrrole structures matches and exclude cyclic porphyrins
    if mol.HasSubstructMatch(polypyrrole_pattern) and not mol.HasSubstructMatch(cyclic_pyrrole_exclusion):
        return True, "Contains polymeric chain of pyrrole units (linked via carbon atoms)"

    # If pyrrole present but not in polypyrrole form
    if mol.HasSubstructMatch(Chem.MolFromSmarts('n1cccc1')):
        return False, "Contains pyrrole ring, but not in form of a polypyrrole chain"

    return False, "No polypyrrole-defining structure found: lacks polymeric pyrrole chain"