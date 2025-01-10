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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # 1. Ester linkage pattern (for acyl chains at sn-1 and sn-2 positions)
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C;H2]')
    
    # 2. Glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts('O[C;H2][C;H][C;H2]O')
    
    # 3. Phosphate group connected to L-serine
    phosphate_serine_pattern = Chem.MolFromSmarts('P(=O)(O)OCC(N)C(=O)O')
    
    # 4. Complete phosphatidylserine pattern
    phosphatidylserine_pattern = Chem.MolFromSmarts('O[C;H2][C;H](O[P](=O)(O)OCC(N)C(=O)O)[C;H2]OC(=O)C')
    
    # Check for ester linkages (acyl chains)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester group(s), expected at least 2 acyl chains at sn-1 and sn-2 positions"
    
    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"
    
    # Check for phosphate group connected to L-serine
    if not mol.HasSubstructMatch(phosphate_serine_pattern):
        return False, "Phosphate group connected to L-serine not found"
    
    # Check for overall phosphatidylserine structure
    if not mol.HasSubstructMatch(phosphatidylserine_pattern):
        return False, "Overall phosphatidylserine structure not found"
    
    return True, "Molecule is a 3-sn-phosphatidyl-L-serine"