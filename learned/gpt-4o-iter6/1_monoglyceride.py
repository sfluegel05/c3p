"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to match the core structure of 1-monoglyceride:
    # The molecule should include a glycerol backbone with an ester linkage at the 1-position
    # This pattern targets an ester group attached to a glycerol, correcting for cheirolectivity
    glycerol_pattern = Chem.MolFromSmarts('C(OCC(=O)[C])CO')  # Less specific, more general pattern
    
    if mol.HasSubstructMatch(glycerol_pattern):
        return True, "Contains a glycerol backbone with acyl linkage at position 1"

    return False, "Does not have a glycerol backbone with acyl linkage at position 1"