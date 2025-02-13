"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with stereo center (R configuration implies [C@H])
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO)OC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone or incorrect configuration found"
    
    # Check for presence of a phosphate group with ethanolamine
    phosphate_ethanolamine_pattern = Chem.MolFromSmarts("COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(phosphate_ethanolamine_pattern):
        return False, "No phosphate group with ethanolamine found"
    
    # Check for 1-O-acyl group (ester linkage)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No 1-O-acyl ester linkage found"
    
    # Assuming we passed all checks
    return True, "Molecule fits all structural criteria for 1-acyl-sn-glycero-3-phosphoethanolamine"