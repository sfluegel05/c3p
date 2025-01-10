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

    # Look for glycerol backbone with phosphate group pattern
    sn_glycerol_3_phosphate_pattern = Chem.MolFromSmarts("O[C@H]([C@@H](O)COP(=O)(O)O)")  # glycerol 3-phosphate
    if not mol.HasSubstructMatch(sn_glycerol_3_phosphate_pattern):
        return False, "No glycerol 3-phosphate backbone found"
    
    # Look for a fatty acid ester in connection with the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("OC(=O)")  # Represents ester linkage
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    esterified = False
    # Confirm an ester linkage is part of the monoacyl-sn-glycerol structure
    for match in ester_matches:
        ester_bond = mol.GetBondBetweenAtoms(match[0], match[1])
        if ester_bond and not ester_bond.IsInRing():  # Avoiding cyclic esters
            # Determine if this ester linkage involves the glycerol backbone
            if any(mol.HasSubstructMatch(sn_glycerol_3_phosphate_pattern)):
                esterified = True
                break

    # Check for exactly one acyl chain as determined by a single ester linkage
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, expected exactly 1 attached to backbone"

    if not esterified:
        return False, "Acyl chain via ester linkage is not correctly bound to the glycerol backbone"
    
    return True, "Contains glycerol 3-phosphate backbone and one acyl group esterified to it"