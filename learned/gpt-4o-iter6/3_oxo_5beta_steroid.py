"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is any 3-oxo steroid that has a beta-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a simpler SMARTS pattern for the 3-oxo group within a steroid backbone
    oxo_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found in the molecule"
    
    # Define a general pattern for steroid backbones, ensuring rings and connections
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCCC4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core not found"

    # Check for 5beta stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    # Look specifically for the stereocenter at position 5 which should be beta ('R' in most cases)
    stereochemically_correct = any(idx == 5 and code == 'R' for idx, code in chiral_centers)
    
    if not stereochemically_correct:
        return False, "5beta stereochemistry not resolved"

    return True, "Molecule is identified as a 3-oxo-5beta-steroid with appropriate stereochemistry"