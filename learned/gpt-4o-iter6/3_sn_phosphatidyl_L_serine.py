"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule corresponds to 3-sn-phosphatidyl-L-serine based on its SMILES string.
    This class is defined by having acyl substituents at sn-1 and sn-2 positions and a phosphoserine group
    at the sn-3 position of a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refine glycerol sn-3 pattern to improve stereo precision
    glycerol_pattern = Chem.MolFromSmarts("O[C@H](COP(O)(=O)O)[C@@H](OC(=O)C)C")  # sn configuration

    # Match glycerol with precise stereo and acyl recognition on sn-1 and sn-2
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No proper glycerol backbone with sn-3 stereochemistry found"

    # Check for acyl substituents
    acyl_pattern = Chem.MolFromSmarts("C(=O)OC")  # Recognize ester linked fatty acid
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 2"

    # Check for phosphoserine at sn-3 position
    phosphoserine_pattern = Chem.MolFromSmarts("COP(O)(=O)O[C@@H](N)C(=O)O")  # sn-stereochemistry included
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # If all patterns match, classify as 3-sn-phosphatidyl-L-serine
    return True, "Contains the structures indicative of a 3-sn-phosphatidyl-L-serine"