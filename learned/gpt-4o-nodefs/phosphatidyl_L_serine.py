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

    # SMARTS pattern for glycerol backbone with 2 ester bonds
    glycerol_pattern = Chem.MolFromSmarts("C(COC(=O))C(O)COC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two ester bonds found"

    # SMARTS pattern for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # SMARTS pattern for L-serine head group (amino acid moiety)
    serine_pattern = Chem.MolFromSmarts("NCC(C(O)=O)OP")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No L-serine head group found attached to phosphate"

    return True, "Contains phosphatidyl-L-serine structure"