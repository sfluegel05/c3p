"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A molecule is classified as such if it contains a glycerol backbone with specific stereochemistry,
    acyl substituents at the sn-1 and sn-2 positions, and a phosphoserine at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone with specified stereochemistry (sn-3)
    glycerol_sn3_pattern = Chem.MolFromSmarts("OC[C@@H](O[P](=O)(O)OC[C@H](N)C(=O)O)[C@H](OC(=O))C(=O)")
    if not mol.HasSubstructMatch(glycerol_sn3_pattern):
        return False, "No proper glycerol sn-3 stereochemistry found"

    # Look for acyl groups on sn-1 and sn-2
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 2:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 2"

    # Look for phosphoserine moiety
    phosphoserine_pattern = Chem.MolFromSmarts("COP(=O)(OCC(N)C(=O)O)O")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    return True, "Contains the structures indicative of a 3-sn-phosphatidyl-L-serine"