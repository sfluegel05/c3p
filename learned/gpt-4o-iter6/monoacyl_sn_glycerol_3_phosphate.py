"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Glycerol backbone pattern: [O][C@@H](O)COP
    glycerol_backbone = Chem.MolFromSmarts("[O][C@@H](O)COP")
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone with phosphate found"
        
    # Acyl group pattern: [O][C](=O)[C,C][C,C] (attached to glycerol backbone)
    acyl_pattern = Chem.MolFromSmarts("[O][C](=O)[C,C][C,C]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached"

    # Check for a single acyl group
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 1"
    
    return True, "Contains a glycerol backbone with phosphate and one acyl group attached correctly"