"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    This includes a glycerol backbone with an acyl chain at the 1-hydroxy position and a phosphoserine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2][CX4][OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphoserine pattern
    phosphoserine_pattern = Chem.MolFromSmarts("P(=O)(O)O[C@H](N)C(O)=O") 
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # Look for acyl group pattern linked to glycerol
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@H](O)C")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl chain linked to the glycerol 1-hydroxy position"

    return True, "Contains glycerol backbone with acyl chain and attached phosphoserine group"