"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: CHEBI:51601 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.EState import EState_VSA
from rdkit.Chem.GroupContribution import GCVPC

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound based on its SMILES string.
    A volatile organic compound is any organic compound having an initial boiling point
    less than or equal to 250 °C (482 °F) measured at a standard atmospheric pressure
    of 101.3 kPa.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a volatile organic compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for organic compound (contains carbon atoms)
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) == 0:
        return False, "No carbon atoms found, not an organic compound"
    
    # Calculate boiling point using group contribution method
    boiling_point = GCVPC.getGroupContrib(mol)
    
    # Check if boiling point is <= 250 °C
    if boiling_point <= 250:
        return True, f"Estimated boiling point of {boiling_point:.2f} °C is <= 250 °C"
    else:
        return False, f"Estimated boiling point of {boiling_point:.2f} °C is > 250 °C"