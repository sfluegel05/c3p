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
    
    # Identify aldoxime group (C=N-OH) linked to aliphatic carbon
    aldoxime_pattern = Chem.MolFromSmarts('[CH3,CH2,CH,CH0][CX3]=NO')
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No aliphatic aldoxime group found"
    
    # Check for absence of aromatic bonds around the oxime carbon
    aromatic_context_pattern = Chem.MolFromSmarts('[c][CX3](=NO)')
    if mol.HasSubstructMatch(aromatic_context_pattern):
        return False, "Aromatic context detected around oxime, not strictly aliphatic"

    return True, "Contains aliphatic structure with aldoxime group"

# Example usage:
# result, reason = is_aliphatic_aldoxime('C(CCCCCCCSC)=NO')
# print(result, reason)