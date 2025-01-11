"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime has the -C=N-OH group and is derived from non-aromatic (aliphatic) aldehydes.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify the aldoxime group pattern (C=N-OH)
    aldoxime_pattern = Chem.MolFromSmarts('[CX3]=N[OH]')
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No aldoxime group found"
    
    # Ensure the aldoxime group is part of an aliphatic aldehyde
    aliphatic_aldoxime_pattern = Chem.MolFromSmarts('[CX3;!r]=N[OH]')
    if not mol.HasSubstructMatch(aliphatic_aldoxime_pattern):
        return False, "No aliphatic aldehyde origin for aldoxime"
    
    # Check for any aromatic atoms or rings anywhere in the molecule
    if mol.HasSubstructMatch(Chem.MolFromSmarts('a')):
        return False, "Aromatic structures detected in the molecule"

    # Additional check: Ensure the adjacent to the C=N-OH group is aliphatic
    adj_aliphatic_check = Chem.MolFromSmarts('[CX3]([!r])=N[OH]')
    if not mol.HasSubstructMatch(adj_aliphatic_check):
        return False, "Neighboring structures to C=N-OH are not entirely aliphatic"
    
    return True, "Contains aliphatic structure with aldoxime group"

# Example usage:
# result, reason = is_aliphatic_aldoxime('C(CCCCCCCSC)=NO')
# print(result, reason)