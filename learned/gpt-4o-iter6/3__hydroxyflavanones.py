"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone has a flavanone structure with a hydroxy group at the 3' position of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the general flavanone structure with SMILES or SMARTS
    # Basic flavanone pattern (benzopyranone core plus attached phenyl)
    flavanone_pattern = Chem.MolFromSmarts("O=C1[C@@H]([C@H](O)C2=CC=CC=C2)CC3=CC=CC=C3O1")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Flavanone core structure not found"

    # Check for hydroxy substitution at the 3' position of the phenyl group
    hydroxy_3_prime_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    if not mol.HasSubstructMatch(hydroxy_3_prime_pattern):
        return False, "No hydroxy group at 3' position of phenyl ring"

    return True, "Contains flavanone structure with a hydroxy group at 3' position of the phenyl ring"