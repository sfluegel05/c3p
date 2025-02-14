"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoserine has a glycerol backbone with an acyl group at the 1-position,
    and a phosphoserine group at the 3-position.

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

    # Check for glycerol backbone with sn-1 and sn-3 stereochemistry
    # Specifically looking for stereochem. of -C[C@H](O)C- 
    glycerol_pattern = Chem.MolFromSmarts("[CX4][C@H](O)[CX4]")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
        return False, "No sn-glycerol backbone found"
    
    glycerol_c1_idx = glycerol_match[0] # Get the index of the first carbon of glycerol
    glycerol_c2_idx = glycerol_match[1] # Get the index of the center carbon of glycerol
    glycerol_c3_idx = glycerol_match[2] # Get the index of the third carbon of glycerol

    # Check for ester group attached to sn-1 carbon
    # Looking for C-O-C=O
    ester_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    found_sn1_ester = False
    for match in ester_matches:
        # Get the atom connected to the glycerol_c1. It should be the oxygen of the ester
      
      for atom in mol.GetAtomNeighbors(glycerol_c1_idx):
        if atom.GetIdx() == match[1]:
            found_sn1_ester = True
            break
      if found_sn1_ester:
        break
    if not found_sn1_ester:
        return False, "No ester group at sn-1 position"


    # check if chain connected to ester is long enough (>=4 carbons)
    acyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    
    found_acyl_chain = False;

    for acyl_match in acyl_matches:
      for match in ester_matches:
        if mol.GetBondBetweenAtoms(acyl_match[0], match[2]): # Check if first carbon is connected to carbonyl
          found_acyl_chain = True
          break
      if found_acyl_chain:
        break
    if not found_acyl_chain:
      return False, "No acyl chain connected to the ester at sn-1 position or chain is too short"

    # Check for phosphoserine group at sn-3 carbon
    phosphoserine_pattern = Chem.MolFromSmarts("P(=O)(O)-O[CX4][CX4](N)C(=O)O")
    phosphoserine_matches = mol.GetSubstructMatches(phosphoserine_pattern);
    found_ps = False;
    for match in phosphoserine_matches:
        for atom in mol.GetAtomNeighbors(glycerol_c3_idx):
            if atom.GetIdx() == match[2]: # check that the O attached to the P is attached to the glycerol C3
                found_ps = True;
                break
        if found_ps:
            break

    if not found_ps:
      return False, "No phosphoserine group found at sn-3 position"


    return True, "Matches 1-acyl-sn-glycero-3-phosphoserine criteria"