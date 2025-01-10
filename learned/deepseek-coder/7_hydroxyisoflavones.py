"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: CHEBI:77064 7-hydroxyisoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone is an isoflavone with a hydroxy group at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the isoflavone core pattern with a hydroxy group at the 7-position
    # The pattern matches the benzopyran-4-one structure with a phenyl group at the 2-position
    # and a hydroxy group at the 7-position.
    isoflavone_pattern = Chem.MolFromSmarts("[O]=[C]1[C]=[C][C]([OH])=[C][C]=[C]1-c2ccccc2")
    
    # Check if the molecule matches the isoflavone core pattern with a hydroxy group at the 7-position
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Molecule does not match the 7-hydroxyisoflavone core structure"
    
    return True, "Molecule contains the 7-hydroxyisoflavone core structure with a hydroxy group at the 7-position"