"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Generic flavanone structure SMARTS pattern 
    flavanone_pattern = Chem.MolFromSmarts("O=C1C=Cc2ccccc2O[C@@H]1c1ccccc1")

    # Check if the molecule has the flavanone structure
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone backbone present"

    # Specific 4'-hydroxy pattern
    four_prime_hydroxy_pattern = Chem.MolFromSmarts("c1(O)cccc1")
    
    # Find all phenyl rings and check for -OH at the 4' position
    ring_matches = mol.GetSubstructMatches(four_prime_hydroxy_pattern)
    if not ring_matches:
        return False, "4'-hydroxy group not found"
    
    return True, "Contains flavanone backbone with a 4'-hydroxy group"