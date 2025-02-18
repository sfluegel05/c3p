"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone has a hydroxy group at the 4' position on the B-ring of a flavanone core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Improved SMARTS pattern for the flavanone core: benzopyranone (Benzopyran-4-one)
    flavanone_core = Chem.MolFromSmarts("O=C1CC2C(C=CC=C2)Oc3ccccc13")
    
    # SMARTS pattern for 4'-hydroxy group attached to the phenyl ring
    hydroxy_position = Chem.MolFromSmarts("Oc1ccc(cc1)-C2Oc3ccccc3C2=O")

    # Check for flavanone core
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "No flavanone core structure found"

    # Check for hydroxy group at the 4' position
    if not mol.HasSubstructMatch(hydroxy_position):
        return False, "Hydroxy group not found at the 4' position"

    return True, "Matches the 4'-hydroxyflavanone structure"