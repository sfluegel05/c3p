"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    
    A monoacyl-sn-glycerol 3-phosphate has a glycerol backbone, a phosphate group on carbon 3,
    and a single acyl group on carbon 1 or 2.

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

    # Revised Glycerol-backbone with phosphate pattern: [O][C@@H](O)COP(=O)(O)O
    # This accounts for the glycerol backbone and phosphate binding
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[O][CH](C(O)COP(=O)(O)O)C")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate correctly configured"

    # Acyl group pattern for positions 1 or 2: [C](=O)[O][CH]
    # Looking for acyl group attached to glycerol at either position 1 or 2
    acyl_pattern = Chem.MolFromSmarts("[C](=O)[O][CH]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 1"
    
    return True, "Contains a glycerol backbone with phosphate and one acyl group attached correctly"

# Example usage:
# result, reason = is_monoacyl_sn_glycerol_3_phosphate("P(OCC(O)COC(=O)CCCCCCCCC/C=C\C/C=C\CCCCC)(O)(O)=O")
# print(result, reason)