"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
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
    
    # Flavone core pattern: 2-aryl-1-benzopyran-4-one
    flavone_pattern = Chem.MolFromSmarts("c1cc(c2c(o1)c(=O)cc(c2)-c3ccccc3)")
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains 2-aryl-1-benzopyran-4-one skeleton"
    else:
        return False, "Does not contain flavone skeleton"