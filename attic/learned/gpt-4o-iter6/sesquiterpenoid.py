"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdchem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are characterized by a C15 backbone that may be rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # A sesquiterpenoid has 15 carbons, which may have rearrangements or modifications
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 18:  # Allow range to accommodate rearrangements
        return False, f"Expected around 15 carbons, got {c_count}"

    # Check for common sesquiterpenoid substructures or features 
    # (such as a diverse set of cyclic and acyclic arrangements)
    patterns = [
        Chem.MolFromSmarts('C1CCC(C)C1'),  # Simple ring structure
        Chem.MolFromSmarts('C=O'),  # Carbonyl groups
        Chem.MolFromSmarts('C=C'),  # Double bonds indicating unsaturation
        Chem.MolFromSmarts('[OH]'),  # Hydroxy group
        Chem.MolFromSmarts('O=C(O)'),  # Ester groups
    ]

    matches = []
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            matches.append(pattern)

    if len(matches) >= 2:  # Consider it a sesquiterpenoid if it has multiple matching patterns
        return True, "Contains multiple sesquiterpenoid characteristic features"

    return False, "Does not exhibit enough sesquiterpenoid characteristic features"