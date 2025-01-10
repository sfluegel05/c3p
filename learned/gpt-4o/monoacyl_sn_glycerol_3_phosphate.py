"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

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

    # Look for glycerol backbone with phosphate group
    sn_glycerol_3_phosphate_pattern = Chem.MolFromSmarts("C(COP(=O)(O)O)O")  # Part of glycerol + phosphate
    if not mol.HasSubstructMatch(sn_glycerol_3_phosphate_pattern):
        return False, "No glycerol 3-phosphate backbone found"
    
    # Look for a single acyl chain bound as ester
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")  # Represents ester linkage
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl linkages, need exactly 1"

    # Confirming its binding to the glycerol backbone but excluding second acyl group
    esterified = False
    for match in acyl_matches:
        atom_indices = [match[0], match[1]]
        acyl_bonded_to_glycerol = any(mol.GetBondBetweenAtoms(i, j).IsInRing() == False for i in atom_indices for j in mol.GetSubstructMatches(sn_glycerol_3_phosphate_pattern))
        if acyl_bonded_to_glycerol:
            esterified = True
            break

    if not esterified:
        return False, "Acyl chain is not correctly bound to the glycerol backbone via ester linkage"
    
    return True, "Contains glycerol 3-phosphate backbone and one acyl group esterified to it"