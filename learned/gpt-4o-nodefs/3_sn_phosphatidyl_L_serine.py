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
    
    # Specific glycerol backbone with stereochemistry for 3-sn-glycerol
    glycerol_stereo_pattern = Chem.MolFromSmarts("O[C@@H](COP(O)(=O)O[C@H](N)C(=O)O)COC(=O)")
    if not mol.HasSubstructMatch(glycerol_stereo_pattern):
        return False, "Missing specific 3-sn-glycerol backbone with correct stereochemistry"
    
    # Specific phosphate linkage to serine
    phosphate_serine_stereo_pattern = Chem.MolFromSmarts("O[P](=O)(O)OC[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(phosphate_serine_stereo_pattern):
        return False, "No specific phosphate-serine linkage found with the correct stereochemistry"
    
    # Check for ester groups with specific stereochemistry
    ester_pattern = Chem.MolFromSmarts("O[C@H](C)C(=O)[C@H](O)COC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages with the correct stereochemistry, need at least 2"
    
    return True, "Contains 3-sn-phosphatidyl-L-serine structure with specific glycerol, phosphate-serine, and fatty acids"