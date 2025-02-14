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
    sn1_ester_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX3](=[OX1])") 
    sn1_ester_matches = mol.GetSubstructMatches(sn1_ester_pattern)
    
    found_sn1_ester = False;
    for match in sn1_ester_matches:
        #Check that one end of the ester is attached to glycerol C1
      if mol.GetBondBetweenAtoms(glycerol_c1_idx,match[0]):
        found_sn1_ester = True
        sn1_ester_c1_idx = match[0]
        break
    if not found_sn1_ester:
         return False, "No ester group at sn-1 position"

    # check if chain connected to ester is long enough (>=4 carbons)
    acyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    found_acyl_chain = False;
    for atom in mol.GetAtomWithIdx(sn1_ester_c1_idx).GetNeighbors():
      if atom.GetIdx() != glycerol_c1_idx:
         if mol.GetSubstructMatch(acyl_chain_pattern,atom.GetIdx()):
              found_acyl_chain = True;
              break;
    if not found_acyl_chain:
        return False, "No acyl chain connected to the ester at sn-1 position or chain is too short"


    # Check for phosphoserine group at sn-3 carbon
    phosphoserine_pattern = Chem.MolFromSmarts("P(=O)(O)-O[CX4][CX4](N)C(=O)O")
    phosphoserine_matches = mol.GetSubstructMatches(phosphoserine_pattern);
    found_ps = False;
    for match in phosphoserine_matches:
        # Check if phosphate is attached to glycerol C3
        if mol.GetBondBetweenAtoms(glycerol_c3_idx, match[2]):
            found_ps = True;
            break;
    if not found_ps:
      return False, "No phosphoserine group found at sn-3 position"


    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if n_count != 1:
         return False, f"Must have exactly 1 nitrogen, got {n_count}"

    if p_count != 1:
         return False, f"Must have exactly 1 phosphorus, got {p_count}"
    if o_count < 7:
        return False, f"Must have at least 7 oxygens for 3-phosphoserine, got {o_count}"

    if c_count < 6:
        return False, f"Too few carbons, got {c_count}"

    return True, "Matches 1-acyl-sn-glycero-3-phosphoserine criteria"