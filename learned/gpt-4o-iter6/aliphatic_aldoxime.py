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
    
    # Identify aliphatic aldoxime group (C=N-OH) linked to aliphatic carbon
    # Ensure the context is aliphatic by excluding any aromatic atoms
    aldoxime_pattern = Chem.MolFromSmarts('[C;!r][CX3]=N[OH]')
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No aliphatic aldoxime group found"
    
    # Verify absence of aromatic bonds directly involving the oxime group
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[c][CX3](=N[OH])')):
        return False, "Aromatic context detected around oxime, not strictly aliphatic"

    return True, "Contains aliphatic structure with aldoxime group"

# Example usage:
# result, reason = is_aliphatic_aldoxime('C(CCCCCCCSC)=NO')
# print(result, reason)