"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone possesses the isoflavone skeleton with a hydroxy group specifically at the 7-position.

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

    # Define key features of isoflavones (2-phenyl-benzopyran-4-one structure)
    phenyl_ring_pattern = Chem.MolFromSmarts("c1ccccc1")  # Phenyl ring
    benzopyran_pattern = Chem.MolFromSmarts("O=c1cc2ccccc2oc1")  # Core benzopyran-4-one

    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No benzopyran-4-one core found"

    if not mol.HasSubstructMatch(phenyl_ring_pattern):
        return False, "No attached phenyl ring found"

    # Check for hydroxy group presence anywhere (as a soft initial filter)
    hydroxy_pattern = Chem.MolFromSmarts("O[cR2]")  # Generic hydroxy pattern on an aromatic carbon
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    if not hydroxy_matches:
        return False, "No hydroxy group found"

    # Further, verify if one hydroxy is in the 7-position relative to the benzopyran
    seven_hydroxy_pattern = Chem.MolFromSmarts("Oc1cc2ccco2c(=O)c1")
    if not mol.HasSubstructMatch(seven_hydroxy_pattern):
        return False, "Hydroxy group not found at the 7-position relative to benzopyran"

    return True, "Has isoflavone skeleton with hydroxy group at 7-position"