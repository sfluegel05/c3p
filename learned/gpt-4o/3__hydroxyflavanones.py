"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent at position 3' of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core pattern (this is a general struct for flavanones)
    flavanone_pattern = Chem.MolFromSmarts("O=C1C(OC2=CC=CC=C2)C[C@H](O)C1")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone skeleton found"

    # Define SMARTS pattern for the hydroxy group at the 3' position of the B-ring
    three_prime_hydroxy_pattern = Chem.MolFromSmarts("Oc1cc([A])ccc1")

    if not mol.HasSubstructMatch(three_prime_hydroxy_pattern):
        return False, "No hydroxy group found exactly at position 3' on the B-ring"

    return True, "Contains flavanone skeleton with a hydroxy group at the 3' position on the phenyl ring"