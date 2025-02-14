"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoserine is an sn-glycerophosphoserine compound having 
    an acyl substituent at the 1-hydroxy position.

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

    # Define SMARTS pattern for the glycerol backbone
    # We remove chiral specifications to be more inclusive
    glycerol_smarts = "[O][C](CO)(CO)"  # Central carbon with two hydroxymethyl groups

    # Define SMARTS pattern for the acyl group at sn-1 position (ester linkage)
    acyl_sn1_smarts = "[C](OC(=O)[C])"  # Carbon with ester linkage

    # Define SMARTS pattern for the phosphoserine group at sn-3 position
    phosphoserine_smarts = "[O]P(=O)([O])[O][C](CO)[N]C(=O)O"  # Phosphoserine moiety

    # Convert SMARTS to RDKit molecule objects
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    acyl_sn1_pattern = Chem.MolFromSmarts(acyl_sn1_smarts)
    phosphoserine_pattern = Chem.MolFromSmarts(phosphoserine_smarts)

    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"

    # Check for acyl group at sn-1 position
    if not mol.HasSubstructMatch(acyl_sn1_pattern):
        return False, "Acyl group not found at sn-1 position"

    # Check for phosphoserine group at sn-3 position
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Phosphoserine group not found at sn-3 position"

    return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoserine"