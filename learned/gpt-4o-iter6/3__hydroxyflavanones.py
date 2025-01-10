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

    # Define the general flavanone core structure
    # Ar-C(=O)C1C(O)C=CC1 with a specific ring structure.
    flavanone_pattern = Chem.MolFromSmarts("O=C1CC2=CC=CC(O)=C2C1C3=CC=CC(O)=C3")
    # Capable of identifying different derivations and stereochemistry forms
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Flavanone core structure not found"

    # Check for hydroxy substitution at the 3' position of the phenyl group
    # Correcting substitution pattern: hydroxy at the meta position relative to C-connection
    hydroxy_3_prime_pattern = Chem.MolFromSmarts("c1(c(O)ccc1)cc(O)c(O)c1")
    if not mol.HasSubstructMatch(hydroxy_3_prime_pattern):
        return False, "No hydroxy group at 3' position of phenyl ring"

    return True, "Contains flavanone structure with a hydroxy group at 3' position of the phenyl ring"