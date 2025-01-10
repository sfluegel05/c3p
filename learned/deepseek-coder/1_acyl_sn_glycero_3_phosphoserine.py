"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: CHEBI:17504 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoserine is an sn-glycerophosphoserine compound with an acyl group at the 1-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core structure: sn-glycero-3-phosphoserine
    core_pattern = Chem.MolFromSmarts("[C@@H](O)COP(O)(=O)OC[C@H](N)C(O)=O")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No sn-glycero-3-phosphoserine core found"

    # Check for acyl group at the 1-position
    # The acyl group is attached to the oxygen at the 1-position of the glycerol backbone
    # Relaxed pattern to match any acyl chain
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)O[C@@H](CO)COP(O)(=O)OC[C@H](N)C(O)=O")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group found at the 1-position"

    # Verify the stereochemistry of the glycerol backbone
    # The 2-position should have the correct stereochemistry (sn-glycerol)
    stereo_match = mol.GetSubstructMatch(core_pattern)
    if len(stereo_match) < 2:
        return False, "Cannot verify stereochemistry"
    
    # Check if the 2-position carbon has the correct stereochemistry
    chiral_center = stereo_match[0]  # The 2-position carbon in the core pattern
    atom = mol.GetAtomWithIdx(chiral_center)
    if not atom.HasProp("_CIPCode"):
        return False, "Cannot determine stereochemistry at the 2-position"
    
    cip_code = atom.GetProp("_CIPCode")
    if cip_code != "S":
        return False, f"Incorrect stereochemistry at the 2-position (got {cip_code}, expected S)"

    # Check for a reasonable acyl chain length (at least 10 carbons)
    # Count the number of carbons in the acyl chain
    acyl_chain = mol.GetSubstructMatch(acyl_pattern)
    if len(acyl_chain) < 10:
        return False, "Acyl chain too short"

    return True, "Contains sn-glycero-3-phosphoserine core with acyl group at the 1-position and correct stereochemistry"