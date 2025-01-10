"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    This compound has a glycerol backbone with one acyl group at positions 1 or 2, and a phosphate at position 3.

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

    # Look for glycerol backbone pattern with phosphate group
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate group found"
    
    # Look for a single acyl chain
    # Note: We assume the acyl chain is a long carbon chain esterified to the glycerol at position 1 or 2
    acyl_pattern = Chem.MolFromSmarts("O=C(O)C")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl chains, need exactly 1"

    # Check for correct stereochemistry at the glycerol's second carbon
    stereo_center = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.CHI_UNSPECIFIED]
    if not stereo_center:
        return False, "Stereochemistry not specified for the glycerol backbone"

    return True, "Contains glycerol backbone with one acyl group and a phosphate group"