"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes, which are C10 skeletons,
    potentially with rearrangements or small modifications like the loss of skeletal atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general C10-terpene like structure (2 isoprenoids pattern)
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[CH][CH2][CH3]")

    # Complex rearrangements and skeleton checking
    terpene_like_pattern = Chem.MolFromSmarts("CC(C)C=CC=CCC(C)C")
    
    if mol.HasSubstructMatch(isoprene_pattern) or mol.HasSubstructMatch(terpene_like_pattern):
        return True, "Structure matches typical patterns seen in monoterpenoids"
    
    # Attempt checking C10 number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count >= 9 and c_count <= 11: 
        # Typically, monoterpenoids have rearranged C10 but can vary from 9 to 11 due to modifications
        return True, "Structure has monoterpene backbone or close to C10 arrangement"

    return False, "No typical monoterpenoid structural features detected"