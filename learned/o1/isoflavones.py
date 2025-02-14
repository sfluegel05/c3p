"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: Isoflavones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-1-benzopyran-4-one skeleton (3-aryl-4H-chromen-4-one)
    and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to handle aromaticity properly
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Chem.KekulizeException:
        return False, "Kekulization failed"

    # Define SMARTS pattern for the isoflavone core
    # Isoflavone core: 3-phenylchromen-4-one
    isoflavone_pattern = Chem.MolFromSmarts("""
    [$([#6]-1=[$([#6](-[#6])=[#6]-[#6]=[#6]-[#8]-1-[#6]=O)])]-c1ccccc1
    """)

    if isoflavone_pattern is None:
        return False, "Error in SMARTS pattern"

    # Check if the molecule matches the isoflavone core pattern
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Does not contain the 3-aryl-1-benzopyran-4-one skeleton"

    # Check for substitutions on the core (which are allowed)
    # So we can accept the match as isoflavone
    return True, "Contains the 3-aryl-1-benzopyran-4-one skeleton characteristic of isoflavones"