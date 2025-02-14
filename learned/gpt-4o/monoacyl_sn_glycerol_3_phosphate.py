"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.

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

    # Look for glycerol backbone with phosphate substitution
    glycerol_pattern = Chem.MolFromSmarts("OC[C@H](O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with phosphate group found"
    
    # Look for a single ester group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check if there is a long carbon chain characteristic of fatty acids
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)O[C@H](O)")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) != 1:
        return False, "No single lipophilic chain typical of an acyl group found"

    # Ensure there is exactly one glycerol phosphate backbone
    phosphate_match = [match for match in fatty_acid_matches if mol.GetAtomWithIdx(match[1]).GetDegree() == 3]
    if len(phosphate_match) != 1:
        return False, "Not exactly one glycerol phosphate backbone"

    return True, "Contains glycerol backbone with one acyl group and a phosphate group"