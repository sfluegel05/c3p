"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavonoid, specifically a flavone, 
    based on its SMILES string. 
    A flavone is characterized by a 2-aryl-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generalized flavone core pattern, allowing substitutions on the rings
    flavone_pattern = Chem.MolFromSmarts("c1cc(C=c2oc3cc(ccc3oc2=O)c2ccccc2)ccc1")
    # Note: We're accounting for substitution flexibility on relevant positions.
    
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains generalized 2-aryl-1-benzopyran-4-one skeleton"
    else:
        return False, "Does not contain flavone skeleton"