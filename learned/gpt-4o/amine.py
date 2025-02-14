"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is derived from ammonia by replacing one, two, or three hydrogen atoms with hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse Mol from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Common patterns to detect amines
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2][!C=O]")  # A non-carbonyl nitrogen
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1][!C=O,#1][!C=O,#1]")  # Non-carbonyl, bound to two non-hydrogens
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3][!C=O,#1][!C=O,#1][!C=O,#1]")  # Non-carbonyl, bound to three non-hydrogens
    cyclic_amine_pattern = Chem.MolFromSmarts("[NX3]@[#6]")  # Nitrogen within a ring bonded to carbon

    # Check the molecule for primary, secondary, tertiary, or cyclic amine patterns
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Molecule contains a primary amine group"
    elif mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Molecule contains a secondary amine group"
    elif mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Molecule contains a tertiary amine group"
    elif mol.HasSubstructMatch(cyclic_amine_pattern):
        return True, "Molecule contains a cyclic amine group"

    return False, "No amine group detected"