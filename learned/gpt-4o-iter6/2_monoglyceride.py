"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has an acyl group esterified at the second carbon of a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a glycerol backbone pattern with esterified second position
    glycerol_2_mono_pattern = Chem.MolFromSmarts("[CX4](O)(CO)COC(=O)C")
    
    if not mol.HasSubstructMatch(glycerol_2_mono_pattern):
        return False, "No 2-monoglyceride pattern found, missing esterified second position"

    # Count carbons and verify ester group
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Ester group not found at second position"

    # Check for hydroxyl groups on the first and third carbon positions
    hydroxyl_pattern = Chem.MolFromSmarts("C(O)COC(=O)C")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing hydroxyl groups on the first or third carbon position"

    return True, "Contains 2-monoglyceride structure with acyl group esterified at the second position"