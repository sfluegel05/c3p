"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    This compound contains a glycerol backbone with attached fatty acids and a phosphate group
    connected to L-serine.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns to detect the required features
    # Glycerol backbone stereo configuration (more generalized pattern)
    glycerol_pattern = Chem.MolFromSmarts("[C@H](OC)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing or incorrect glycerol backbone stereochemistry"
    
    # Phosphate-serine connection allowing variation in charge and stereochemistry
    phosphate_serine_pattern = Chem.MolFromSmarts("COP(O)(=O)OC[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(phosphate_serine_pattern):
        phosphate_serine_neg_charge_pattern = Chem.MolFromSmarts("COP(=O)([O-])OC[C@H](N)C(=O)O")
        if not mol.HasSubstructMatch(phosphate_serine_neg_charge_pattern):
            return False, "No phosphate-serine link found"
    
    # Look for ester linkages indicating fatty acid chains with flexibility in chain lengths
    fatty_acid_ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(fatty_acid_ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"
    
    return True, "Contains 3-sn-phosphatidyl-L-serine structure with glycerol, phosphate-serine, and fatty acids"