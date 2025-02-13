"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose skeleton is closely related
    to cholestan-3-ol, with possible additional carbon atoms in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sterol core fingerprint (adjusted to match sterol carbon frameworks related to cholestan)
    # A simple re-check on Cholesterol as a baseline
    cholesterol_smarts = "C[C@H](O)C1CC[C@]2(C)C3CCC4C(C)(CCC4C3CCC2C1)C"
    cholesterol_pattern = Chem.MolFromSmarts(cholesterol_smarts)
    if not mol.HasSubstructMatch(cholesterol_pattern):
        return False, "Does not contain a sterol core structure"

    # Validate the presence of the 3-hydroxy group (must be explicitly on 3rd carbon of the backbone)
    # Ensure hydroxyl is attached in cholesterol's descriptive landmark
    hydroxy_pattern = Chem.MolFromSmarts("C[C@H](O)[C@H]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Does not have a specific 3-hydroxy group"

    # Allow for variations in C31 or related extensions but must maintain sterol's signature shape
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts("C[C@H](O)C"))):
        return True, "Recognized as a sterol-related framework with a 3-hydroxy group"

    return False, "Failed to match sterol's unique skeletal and functional group requirements"