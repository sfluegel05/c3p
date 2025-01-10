"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate includes a glycerol backbone with one acyl chain and a phosphate 
    group attached to the third position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a glycerol backbone (3 carbon chain with 3 oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[OX2][CH2][CH](O[H0])CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"

    # Check for exactly one ester group (-O-C(=O)R)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly one ester group, found {len(ester_matches)}"

    # Ensure presence of phosphate group (P=O(O)(O)O) bonded to glycerol
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found at third position"

    return True, "Contains glycerol backbone, single acyl chain, and phosphate group at the third position"