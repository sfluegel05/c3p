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

    # Flavanone structure: benzene ring attached to a chromane core
    flavanone_pattern = Chem.MolFromSmarts("C1=CC=C2C(=C1)C(=O)CC2-O")
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone backbone found"

    # Define the 3'-hydroxy substitution more broadly
    # By looking for the hydroxy group on a phenyl group that is part of a flavanone system
    hydroxy_3prime_pattern = Chem.MolFromSmarts("c1cc(O)cc(c1)-C2=CC(O)=C(O2)")
    
    if not mol.HasSubstructMatch(hydroxy_3prime_pattern):
        return False, "No 3'-hydroxy substitution on phenyl ring found"
    
    return True, "Contains 3'-hydroxyflavanone structure"