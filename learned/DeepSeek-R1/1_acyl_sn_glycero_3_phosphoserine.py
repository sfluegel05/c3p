"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: CHEBI:XXXXX 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    The molecule must have:
    - sn-glycero-3-phosphoserine backbone
    - One acyl group (ester) at the sn-1 position
    - No other acyl groups

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the criteria, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core sn-glycero-3-phosphoserine pattern
    # Core structure: glycerol with phosphate at sn-3 connected to serine, acyl at sn-1
    core_pattern = Chem.MolFromSmarts("[C@H](O)(COP(=O)(O)OC[C@H](N)C(=O)O)COC(=O)")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core sn-glycero-3-phosphoserine structure not found"
    
    # Check for exactly one ester group (acyl at sn-1)
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, expected exactly 1"
    
    # Verify the ester is at sn-1 (part of the core)
    # Get the matched atoms in the core pattern to find the ester's position
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Core structure match failed"
    
    # The ester in the core is the OC(=O) group attached to the first carbon of glycerol
    # Assuming the core SMARTS has the ester as part of the pattern, the presence of the core ensures correct position
    
    # Check molecular weight (optional)
    # Typical examples have higher MW, but may not be necessary
    
    # Check for absence of additional acyl groups on the glycerol
    # The core pattern already includes the sn-1 acyl, so any other esters would have been counted
    
    return True, "Contains sn-glycero-3-phosphoserine with an acyl group at the 1-position"