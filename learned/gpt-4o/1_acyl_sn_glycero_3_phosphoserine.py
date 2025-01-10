"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    This class is defined by having an acyl group at the sn-1 position of a glycerophosphoserine backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define substructure patterns with improved specificity
    glycerol_pattern = Chem.MolFromSmarts("O[C@H](CO)CO")  # Stereospecific (glycerol backbone)
    acyl_sn1_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]([O])COP(=O)(O)OC[C@H](N)C(=O)O") # Acyl at sn-1 with phosphoserine

    # Check for glycerol pattern
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for acyl group position with phosphoserine
    acyl_sn1_matches = mol.GetSubstructMatches(acyl_sn1_pattern)
    if len(acyl_sn1_matches) == 0:
        return False, "Missing acyl group at sn-1 position, or incorrect positioning relative to phosphoserine"

    return True, "Molecule contains 1-acyl-sn-glycero-3-phosphoserine"