"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone is characterized by a flavanone backbone with a hydroxyl group
    attached at the 3' position of the phenyl ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for the classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more comprehensive flavanone pattern, cautiously broadening scope
    flavanone_structure = Chem.MolFromSmarts("C1=CC=C2C(=C1)CC=C(O2)C(=O)C")
    
    if not mol.HasSubstructMatch(flavanone_structure):
        return False, "No flavanone backbone found"

    # Define 3'-hydroxy substitution more flexibly
    # Assume the traditional '3' position on one benzene ring for a hydroxy 
    hydroxy_3prime_pattern = Chem.MolFromSmarts("c1c(O)cc(cc1)O")
    
    if not mol.HasSubstructMatch(hydroxy_3prime_pattern):
        return False, "No 3'-hydroxy substitution on phenyl ring found"
    
    return True, "Contains 3'-hydroxyflavanone structure"