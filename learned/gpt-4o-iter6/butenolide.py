"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide consists of a 2-furanone skeleton and its substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extend the SMARTS pattern to allow substitutions on the 2-furanone core
    # This pattern matches a furanone, allowing substituents at all positions (R groups)
    # Pattern allows the central lactone and basic furan ring with R substituents
    furanone_pattern = Chem.MolFromSmarts("O=C1OC=CC1[R0,R1,R2,R3,R4]")
    
    # Match furanone skeleton with substitutions
    if mol.HasSubstructMatch(furanone_pattern):
        return True, "Contains 2-furanone skeleton with allowed substitutions"
    else:
        return False, "Does not contain 2-furanone skeleton or allowed substitutions"