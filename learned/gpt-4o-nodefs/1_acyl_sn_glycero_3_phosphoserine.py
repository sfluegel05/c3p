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
        bool: True if the molecule is 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for key structural components:
    # - Glycerol with phosphate and acyl chain (stereospecifically)
    acyl_glycerol_phosphate_pattern = Chem.MolFromSmarts("[C@@H](OC(=O)C)[CH2]O[P](=O)(O[CH][CH](C(O)=O)N)[O]")
    
    # Check for the key structural component patterns
    if not mol.HasSubstructMatch(acyl_glycerol_phosphate_pattern):
        return False, "1-acyl-sn-glycero-3-phosphoserine structure not matched"

    return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"