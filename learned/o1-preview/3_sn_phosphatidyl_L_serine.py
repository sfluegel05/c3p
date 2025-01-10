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
    at the sn-1 and sn-2 positions.

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

    # Define SMARTS pattern for the phosphatidyl-L-serine core with stereochemistry
    phosphatidylserine_smarts = """
    [C@@H](CO[P](=O)(O)OC[C@H](N)C(=O)O)(OC(=O)[#6])[C@@H](OC(=O)[#6])O
    """

    phosphatidylserine_pattern = Chem.MolFromSmarts(phosphatidylserine_smarts)
    if phosphatidylserine_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for the phosphatidyl-L-serine core structure with correct stereochemistry
    if not mol.HasSubstructMatch(phosphatidylserine_pattern, useChirality=True):
        return False, "Molecule does not match the phosphatidyl-L-serine core structure with correct stereochemistry"

    # Optionally, check for two acyl chains (R groups) attached via ester linkages
    acyl_ester_pattern = Chem.MolFromSmarts('OC(=O)[C]')
    acyl_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found {len(acyl_matches)} acyl ester group(s), expected 2 at sn-1 and sn-2 positions"

    return True, "Molecule is a 3-sn-phosphatidyl-L-serine"