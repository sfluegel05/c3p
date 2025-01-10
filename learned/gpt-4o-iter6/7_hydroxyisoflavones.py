"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone possesses a hydroxyisoflavone skeleton with a hydroxy group specifically at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Core scaffold for isoflavones: C2-C3 linkage with a benzopyran ring leading to a phenyl
    isoflavone_core_pattern = Chem.MolFromSmarts("O=C1C=Cc2ccccc2Oc2c1cccc2")
    
    if not mol.HasSubstructMatch(isoflavone_core_pattern):
        return False, "No isoflavone core structure found"
    
    # SMARTS pattern to detect hydroxy group position - flexible for the A-ring
    hydroxy_7_position_pattern = Chem.MolFromSmarts("Oc1ccc2cc(o)cc2c1")
    
    if not mol.HasSubstructMatch(hydroxy_7_position_pattern):
        return False, "No hydroxy group at the 7-position found"
    
    return True, "Has isoflavone skeleton with hydroxy group at 7-position"