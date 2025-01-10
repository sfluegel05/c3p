"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone has a flavanone core structure with a hydroxy group
    on the phenyl ring specifically at the 3' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the flavanone core structure
    flavanone_core_pattern = Chem.MolFromSmarts("O=C1[C@H](Oc2ccccc2)CC2=CC=CC=C12")
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core structure not found"

    # Define a SMARTS pattern for the hydroxy group at the 3' position
    # Assuming the 3' hydroxy position is the meta position relative to the ring attachment
    # The phenyl ring is part of the flavanone core, hence might need reevaluation
    three_prime_hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")

    # Find if there's a phenyl group matching this pattern
    match = mol.HasSubstructMatch(three_prime_hydroxy_pattern)
    
    if not match:
        return False, "No hydroxy group at 3' position on the phenyl ring"
    
    return True, "Contains flavanone core with a hydroxy group at 3' position"