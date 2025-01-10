"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define substructure patterns 
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    phosphoserine_pattern = Chem.MolFromSmarts("OC1COP([O-])([O-])=O")  # phosphate with serine
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")  # Simple acyl ester linkage

    # Check for glycerol pattern (may need to adjust based on stereochemistry specifics)
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Check for phosphoserine group pattern
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # Check for acyl ester linkage at the 1-position (simplified check for one acyl linkage)
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Missing acyl group at sn-1 position"

    return True, "Molecule contains sn-glycero-3-phosphoserine with acyl group at sn-1 position"