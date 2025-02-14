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
    # More flexible SMARTS to find glycerol backbone with a phosphate group.
    glycerol_phosphate_pattern = Chem.MolFromSmarts("C[CH](O)COP(=O)(O)(O)")
    glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if not glycerol_phosphate_matches:
      glycerol_phosphate_pattern = Chem.MolFromSmarts("C[CH](O)C[O]P(=O)(O)(O)") #trying alternative because the phosphate group can be represented with single or double bonds to the O.
      glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
      if not glycerol_phosphate_matches:
        glycerol_phosphate_pattern = Chem.MolFromSmarts("O[CH](CO)C[O]P(=O)(O)(O)")
        glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
      if not glycerol_phosphate_matches:
        return False, "No sn-glycerol-3-phosphate core found"

    # Verify the chirality of the central carbon is S
    glycerol_carbon_index = glycerol_phosphate_matches[0][1]
    chiral_center = mol.GetAtomWithIdx(glycerol_carbon_index)
    if chiral_center.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
         return False, "Central carbon is not S chiral"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    #check ester connected to either position 1 or 2. Not the phosphate
    ester_match = ester_matches[0]
    carbonyl_carbon_index = ester_match[1]
    
    glycerol_phosphate_match = glycerol_phosphate_matches[0]
    glycerol_carbons = [glycerol_phosphate_match[0], glycerol_phosphate_match[1]]
    
    ester_connected = False
    for glycerol_carbon in glycerol_carbons:
      for neighbor in mol.GetAtomWithIdx(glycerol_carbon).GetNeighbors():
        if neighbor.GetIdx() == carbonyl_carbon_index:
            ester_connected = True
            break
      if ester_connected:
        break
    if not ester_connected:
       return False, "Ester not connected to position 1 or 2 of glycerol"

    # Check for long carbon chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing acyl chain"
    
    # count carbons and hydrogens to check for long chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if c_count < 6 or h_count < 10:
        return False, "Not a long chain fatty acid"
    
    return True, "Contains a glycerol backbone, a phosphate group at position 3, and a single acyl group at position 1 or 2"