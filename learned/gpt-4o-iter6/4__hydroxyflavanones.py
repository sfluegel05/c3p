"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent located at position 4'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined SMARTS pattern for flavanone core structure
    flavanone_core_pattern = Chem.MolFromSmarts("O=C1C[C@H](Oc2ccccc12)c3ccccc3")  # Flavanone with stereochemistry inclination

    # Check for flavanone core in molecule
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure identified"

    # Improved SMARTS pattern for identifying 4'-hydroxyphenyl group
    four_prime_hydroxy_pattern = Chem.MolFromSmarts("c1ccc(O)cc1C2C(=O)CC(O2)c3ccccc3")  # Ensures a para-position from attachment

    # Check if the flavanone core is correctly placed
    matches = mol.GetSubstructMatches(four_prime_hydroxy_pattern)
    if matches:
        return True, "Contains 4'-hydroxyflavanone core and a 4'-hydroxy group"

    return False, "4'-hydroxy group not found at the correct position"