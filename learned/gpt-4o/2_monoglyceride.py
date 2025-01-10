"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone with an acyl substituent at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern with C2 being connected to an ester
    # The glycerol pattern can be identified as C(CO)CO, and we need to verify that the ester is attached to the middle carbon
    glycerol_pattern = Chem.MolFromSmarts("OC(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone in the expected form"

    # Look specifically for an acyl chain at the second position via ester linkage (O=C-O)
    ester_pattern = Chem.MolFromSmarts("C(OC(=O)[*])CO")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No acyl substituent at the 2-position via ester linkage"

    return True, "Contains glycerol backbone with an acyl substituent at the 2-position (ester linkage)"