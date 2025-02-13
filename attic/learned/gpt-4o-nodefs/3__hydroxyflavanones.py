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

    # Find flavanone backbone pattern
    flavanone_pattern = Chem.MolFromSmarts("O=C1CC2=C(O1)C=CC=C2")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone backbone found"
    
    # Find 3'-hydroxy substitution on a phenyl group
    hydroxy_3prime_pattern = Chem.MolFromSmarts("c1cc(O)cc(c1)O")
    if not mol.HasSubstructMatch(hydroxy_3prime_pattern):
        return False, "No 3'-hydroxy substitution on phenyl ring found"
    
    return True, "Contains 3'-hydroxyflavanone structure"