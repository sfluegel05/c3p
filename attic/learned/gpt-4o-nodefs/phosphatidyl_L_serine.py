"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for glycerol backbone connected to the fatty acids and phosphate
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)(O)O")  # Simplified pattern
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two ester bonds and a phosphate found"

    # Flexible SMARTS pattern for the L-serine group attached via phosphate
    serine_pattern = Chem.MolFromSmarts("N[C@@H](C(O)=O)COP(=O)O")  # Adjusted for potential variations
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No L-serine head group found attached via phosphate linkage"

    return True, "Contains phosphatidyl-L-serine structure"