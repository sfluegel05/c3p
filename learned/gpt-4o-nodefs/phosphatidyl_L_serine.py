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

    # SMARTS pattern for glycerol backbone linked to esters and phosphate
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO[P](=O)(O)O[C@H](COC(=O)[CX4,CX3])OC(=O)[CX4,CX3])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two ester bonds and a phosphate found"

    # SMARTS pattern for serine attached via the phosphate group
    # Ensure it can catch variations in the phosphatidyl-L-serine moiety
    serine_pattern = Chem.MolFromSmarts("N[C@@H](C(O)=O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No L-serine head group found attached via phosphate linkage"

    return True, "Contains phosphatidyl-L-serine structure"