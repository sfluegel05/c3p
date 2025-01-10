"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define individual patterns for key structural components:
    # - Phosphate group connected to glycerol
    phosphate_glycerol_pattern = Chem.MolFromSmarts("P(=O)(O)OC[C@H](O)CO")

    # - Acyl chain connected to glycerol
    acyl_glycerol_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CO)")

    # - Phosphoserine side chain
    phosphoserine_pattern = Chem.MolFromSmarts("OC[C@H](N)C(O)=O")

    # Check each key structural component pattern in the molecule
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
        return False, "Phosphate group linked to glycerol not found"

    if not mol.HasSubstructMatch(acyl_glycerol_pattern):
        return False, "Acyl chain linkage to glycerol missing"

    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Phosphoserine structure not identified"

    return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"