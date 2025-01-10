"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: CHEBI:17504 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoethanolamine has a glycerol backbone with (R)-configuration,
    a single acyl chain at the 1-position, and a phosphoethanolamine group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with (R)-configuration
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O)(COC(=O)[*])COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with (R)-configuration found"

    # Check for single acyl chain at the 1-position
    acyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 1:
        return False, "No acyl chain found at the 1-position"

    # Check for phosphoethanolamine group at the 3-position
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found at the 3-position"

    # Check for at least one ester bond (single acyl chain)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, f"Found {len(ester_matches)} ester bonds, need at least 1"

    # Check molecular weight - typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for 1-acyl-sn-glycero-3-phosphoethanolamine"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for 1-acyl-sn-glycero-3-phosphoethanolamine"
    if o_count < 6:
        return False, "Too few oxygens for 1-acyl-sn-glycero-3-phosphoethanolamine"

    return True, "Contains glycerol backbone with (R)-configuration, single acyl chain at 1-position, and phosphoethanolamine group at 3-position"