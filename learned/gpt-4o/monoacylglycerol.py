"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (O-C-C-C-O)
    glycerol_backbone_pattern = Chem.MolFromSmarts("OCCO")
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "No glycerol backbone found"
    
    # Look for acyl group pattern (O=C-C)
    acyl_group_pattern = Chem.MolFromSmarts("O=C(C)")
    acyl_groups = mol.GetSubstructMatches(acyl_group_pattern)
    if len(acyl_groups) != 1:
        return False, f"Found {len(acyl_groups)} acyl groups, requires exactly 1"
    
    # Ensure the acyl group is connected to glycerol backbone
    ester_bond_pattern = Chem.MolFromSmarts("O=C(C)OCCO")
    if not mol.HasSubstructMatch(ester_bond_pattern):
        return False, "Acyl group not properly ester linked to glycerol backbone"
    
    # Verify molecular parts: 3 carbons for glycerol, at least 2 oxygens, and one acyl chain
    glycerol_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and len(atom.GetNeighbors()) == 3 and any(nei.GetAtomicNum() == 8 for nei in atom.GetNeighbors()))
    if glycerol_carbons < 3:
        return False, "Too few carbons for a complete glycerol backbone"

    glycerol_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and any(nei.GetAtomicNum() == 6 for nei in atom.GetNeighbors()))
    if glycerol_oxygens < 2:
        return False, "Too few oxygens for a complete glycerol structure"
    
    return True, "Contains a glycerol backbone with one acyl group and two hydroxyl groups"