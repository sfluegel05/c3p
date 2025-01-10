"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    Corrinoids contain a corrin nucleus, which consists of four reduced or partly reduced pyrrole rings 
    joined in a macrocycle by three =C- groups and one direct carbon-carbon bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for corrin nucleus: four pyrrole-like rings linked by three =C- bridges and one direct C-C bond
    # This is an approximation to capture the essence of a corrin structure.
    # Individual pyrrole rings: C1=NC=CC=C1
    # Connecting pyrrole-like units: [C,c]=[N,n][C,c]=[C,c][C,c]=[N,n][C,c]=[C,c][C,c]=[N,n][C,c]=[C,c][C,c]=[N,n][C,c]=[C,c]

    # Corrin-like core pattern (approximation using common connectivity pattern):
    corrin_pattern = Chem.MolFromSmarts('[C,c]1=[N,n][C,c]=[C,c][C,c]=[N,n]2[C,c]=[N,n]3[C,c]=[C,c][C,c]4=[N,n]1C2=[C,c]3[C,c]=[N,n]4')

    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Contains a corrin-like nucleus pattern"
    
    return False, "Does not contain a discernible corrin nucleus pattern"