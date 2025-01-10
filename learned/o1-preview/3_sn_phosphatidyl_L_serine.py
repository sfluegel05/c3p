"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:8339 3-sn-phosphatidyl-L-serine
"""

from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine is a glycerophosphoserine compound having acyl substituents
    at the 1- and 2-hydroxy positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    from rdkit import Chem

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # 1. Acyl ester linkage pattern (for acyl chains at sn-1 and sn-2 positions)
    acyl_ester_pattern = Chem.MolFromSmarts('C(=O)O[C]')
    
    # 2. Glycerol backbone pattern (less restrictive)
    glycerol_pattern = Chem.MolFromSmarts('[O][C][C][C][O]')
    
    # 3. Phosphate connected to serine pattern (generalized)
    phosphate_serine_pattern = Chem.MolFromSmarts('P(=O)(O)OCC(N)C(=O)O')
    
    # Check for at least two acyl ester linkages (sn-1 and sn-2 positions)
    acyl_ester_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    if len(acyl_ester_matches) < 2:
        return False, f"Found {len(acyl_ester_matches)} acyl ester group(s), expected at least 2 acyl chains at sn-1 and sn-2 positions"
    
    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"
    
    # Check for phosphate group connected to serine
    if not mol.HasSubstructMatch(phosphate_serine_pattern):
        return False, "Phosphate group connected to serine not found"
    
    return True, "Molecule is a 3-sn-phosphatidyl-L-serine"