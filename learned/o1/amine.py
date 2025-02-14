"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:32988 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two, or three hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for primary, secondary, and tertiary amines
    # Exclude amides, nitro groups, and other nitrogen-containing groups
    amine_pattern = Chem.MolFromSmarts("[NX3;!$(N-C=O);!$(N=O);!$(N-[!#6]);!$(N=[!#6]);!$(N#[!#6]);!$(N+:*)]")
    if amine_pattern is None:
        return False, "Error in SMARTS pattern"

    # Search for amine groups in the molecule
    matches = mol.GetSubstructMatches(amine_pattern)
    if matches:
        return True, "Molecule contains an amine group"
    else:
        return False, "No amine group found"