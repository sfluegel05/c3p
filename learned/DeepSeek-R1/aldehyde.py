"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is a compound with the structure RC(=O)H, where a carbonyl group is bonded to one hydrogen and one R group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Aldehyde SMARTS pattern: [CH1](=O) where the carbon has exactly one H attached
    aldehyde_pattern = Chem.MolFromSmarts("[CH1]=O")
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if len(matches) == 0:
        return False, "No aldehyde group (RC(=O)H) found"
    elif len(matches) > 1:
        return False, f"Multiple aldehyde groups found ({len(matches)})"
    
    # Check if the aldehyde carbon is terminal (only connected to the carbonyl oxygen and R group)
    # The aldehyde carbon should have exactly two bonds: one to O and one to R
    for match in matches:
        atom = mol.GetAtomWithIdx(match[0])
        if atom.GetDegree() != 2:
            return False, "Aldehyde carbon has incorrect bonding (must be bonded to O and R only)"
    
    return True, "Contains an aldehyde group (RC(=O)H)"