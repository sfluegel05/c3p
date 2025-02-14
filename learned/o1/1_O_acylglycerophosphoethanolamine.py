"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    This class of molecules has a glycerol backbone with:
    - An acyl chain attached via an ester bond at position 1
    - A free hydroxyl group at position 2
    - A phosphoethanolamine group attached at position 3
    - No other acyl chains or substitutions on the glycerol backbone

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for 1-O-acylglycerophosphoethanolamine
    pattern = Chem.MolFromSmarts("""
    [C@@H]1([O][P](=O)([O])[O][C][C][N])  # Chiral carbon at position 3 with phosphoethanolamine
    [C@@H](O)[C@@H](O[C](=O)[C])O1        # Ring closure to represent glycerol with acyl group at position 1
    """)

    # If chiral specification causes issues, use a more general pattern without chiral tags
    if pattern is None:
        # Simplified pattern without chiral centers
        pattern = Chem.MolFromSmarts("""
        O[C@@H](COP(O)(O)OCCN)COC(=O)C
        """)

    if pattern is None:
        return False, "Error in SMARTS pattern"

    # Check for substructure match
    if not mol.HasSubstructMatch(pattern):
        return False, "Molecule does not match 1-O-acylglycerophosphoethanolamine pattern"

    # Ensure there is only one acyl chain attached via ester bond
    ester_pattern = Chem.MolFromSmarts('OC(=O)[C]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, expected exactly 1 at position 1"

    # Check for free hydroxyl group at position 2
    free_oh_pattern = Chem.MolFromSmarts('[C@H](O)[C@H](O)CO[P](=O)(O)OCCN')
    if not mol.HasSubstructMatch(free_oh_pattern):
        return False, "No free hydroxyl group at position 2 found"

    # Ensure there are no other acyl chains
    total_ester_groups = len(ester_matches)
    if total_ester_groups > 1:
        return False, "Additional ester groups found on the molecule"

    return True, "Molecule is a 1-O-acylglycerophosphoethanolamine"