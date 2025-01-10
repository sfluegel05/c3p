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
    
    # An improved SMARTS for corrin nucleus: four pyrrole-like rings linked by three methine (=C-) bridges and one C-C bond
    # The =C- bridges are important, ensuring the correct conjugation and macrocycle formation.
    corrin_pattern = Chem.MolFromSmarts("""
        [cH][c][n][c][n][c]1[c][n][c][n][c]([c][n][c][C]2)
        [cH][c]1[c][n][c][n][c][c]2=[C][C]2[C][C][C]=C1
    """)

    # The above SMARTS is an intricate approximation aiming to represent the unique asymmetric nature of the corrin macrocycle.

    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Contains a corrin-like nucleus pattern"
    
    return False, "Does not contain a discernible corrin nucleus pattern"