"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # Define SMARTS pattern to better match the isoflavone core
    # Capture potential substitutions with 'c' representing aromatic carbons 
    # and variations at the typical substitution sites
    isoflavone_pattern = Chem.MolFromSmarts('Oc1cc2c(c1)cc(=O)oc2-c1ccccc1')
    
    if mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Contains 3-aryl-1-benzopyran-4-one core structure with possible substitutions"

    return False, "Does not contain an isoflavone core structure"