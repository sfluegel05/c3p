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

    # Improved SMARTS pattern for flavanone core structure
    flavanone_core_pattern = Chem.MolFromSmarts("c1ccc2c(c1)O[C@@H](C2=O)Cc3ccccc3")

    # Check for flavanone core in molecule
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure identified"

    # Detect 4'-hydroxy group
    # Based on typical numbering in flavanones where 4'-OH is on the phenyl ring opposite the benzopyranone moiety
    four_prime_hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)ccc1[C@H]2COc3ccccc3C2=O")
    
    if not mol.HasSubstructMatch(four_prime_hydroxy_pattern):
        return False, "4'-hydroxy group not found at the correct position"

    return True, "Contains 4'-hydroxyflavanone core and a 4'-hydroxy group"