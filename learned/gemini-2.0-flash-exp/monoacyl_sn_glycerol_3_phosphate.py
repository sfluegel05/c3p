"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate has a glycerol backbone, a phosphate group at the 3-position,
    and a single fatty acid chain (acyl group) at either position 1 or 2 of the glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sn-glycerol-3-phosphate core with correct stereochemistry
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[C@H]([OX2])[CH2X4][OX2]P(=O)([OX1])([OX1])")
    glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if not glycerol_phosphate_matches:
         glycerol_phosphate_pattern = Chem.MolFromSmarts("[OX2][CH2X4][C@H]([OX2])P(=O)([OX1])([OX1])")
         glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
         if not glycerol_phosphate_matches:
            return False, "No sn-glycerol-3-phosphate core found"
    
    #verify the chirality of the central carbon is S
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    found_correct_chirality = False
    for atom_index, chirality in chiral_centers:
         if chirality == 'S':
            glycerol_carbons = [match[0] for match in glycerol_phosphate_matches]
            glycerol_carbons += [match[1] for match in glycerol_phosphate_matches]
            if atom_index in glycerol_carbons:
               found_correct_chirality = True
    if not found_correct_chirality:
        return False, "Central carbon not S chiral"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"


    #check ester connected to either position 1 or 2. Not the phosphate
    ester_match = ester_matches[0]
    ester_oxygen_atom_index = ester_match[0]
    glycerol_phosphate_match = glycerol_phosphate_matches[0]

    glycerol_carbons_and_oxygens = [glycerol_phosphate_match[0],glycerol_phosphate_match[1], glycerol_phosphate_match[2]]
    for i in mol.GetAtomWithIdx(glycerol_carbons_and_oxygens[0]).GetNeighbors():
      if i.GetIdx() != glycerol_carbons_and_oxygens[1] and i.GetIdx() != glycerol_carbons_and_oxygens[2]:
        if i.GetIdx() != ester_oxygen_atom_index:
          return False, "Ester not connected to position 1 or 2 of glycerol"
    
    for i in mol.GetAtomWithIdx(glycerol_carbons_and_oxygens[2]).GetNeighbors():
      if i.GetIdx() != glycerol_carbons_and_oxygens[0] and i.GetIdx() != glycerol_carbons_and_oxygens[1]:
         if i.GetIdx() != ester_oxygen_atom_index:
           return False, "Ester not connected to position 1 or 2 of glycerol"
    
    
    # Check for long carbon chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing acyl chain"

    # count carbons and hydrogens to check for long chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if c_count < 8 or h_count < 10:
        return False, "Not a long chain fatty acid"
   
    return True, "Contains a glycerol backbone, a phosphate group at position 3, and a single acyl group at position 1 or 2"