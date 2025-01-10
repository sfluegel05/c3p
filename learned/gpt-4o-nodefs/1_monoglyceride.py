"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride typically has a single fatty acid chain attached to one
    of the primary hydroxyl groups of glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a generalized SMARTS pattern for a primary glycerol ester
    # This considers the structure where glycerol is esterified at a primary hydroxyl
    glycerol_ester_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)")
    if not glycerol_ester_pattern:
        return False, "Failed to create glycerol ester SMARTS pattern"

    # Match the glycerol ester pattern in the molecule
    if not mol.HasSubstructMatch(glycerol_ester_pattern):
        return False, "No ester linkage found at a primary position of glycerol"

    # Notice that we don't check for an explicit carbon chain directly attached to the ester
    # This keeps the implementation less prone to misidentifications due to complex branching in SMILES

    return True, "Contains glycerol backbone with an ester linkage at one of the primary positions and a potential fatty acid chain"