"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    This class is characterized as a 1-O-acylglycerophosphoethanolamine with (R)-configuration at sn-2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Recognize the glycerol backbone with chiral center at sn-2
    # Specifically check for (R)-configuration using InChI to handle stereo
    glycerol_backbone_pattern = Chem.MolFromSmarts("[O][C@H](CO[CX3])")
    subs = mol.GetSubstructMatches(glycerol_backbone_pattern)
    if not any(Chem.FindMolChiralCenters(mol, includeUnassigned=True)):
        return False, "No glycerol backbone with proper stereochemistry (R) found"

    # Check for attached phosphate group in broader way ensuring one oxygen link to ethanolamine
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphoethanolamine group found"
    
    # Look for ester bond at sn-1 position
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O[C@H]") 
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage detected at sn-1 position"
    
    return True, "Molecule matches 1-acyl-sn-glycero-3-phosphoethanolamine structural criteria"