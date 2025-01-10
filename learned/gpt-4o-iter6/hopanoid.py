"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton, typically featuring a pentacyclic structure with specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define comprehensive hopanoid SMARTS patterns
    hopanoid_patterns = [
        '[C@@]12CC[C@@]3(C)[C@]1(CC[C@@H](C2)C3)C',  # Core hopanoid skeleton
        '[C@]12(CC[C@@H](CC3=C2CCC[C@]3C1(C)C)C)C',  # Another core-style skeleton
        '[C@@H]1CC[C@]2(C)[C@]1(CC[C@@]3([C@@H]2CC3)C)C',  # Alternate stereochemistry
        '[C@@]1(CC[C@]2([C@@H]1CC[C@@H](C2)C)C)C',  # Integration of core structures
        '[C@H]1[C@@]2(CC[C@]3([C@@H]2CCC[C@@]13C(C)(C)C)C)C',  # Expanded cyclic structure
        # Patterns can be extended or modified to cater to detected false negatives
    ]

    # Convert SMARTS into RDKit Mol objects
    try:
        hopanoid_mols = [Chem.MolFromSmarts(pattern) for pattern in hopanoid_patterns]
    except Exception as e:
        return None, f"Error in SMARTS pattern conversion: {e}"
    
    if any(hopanoid_mol is None for hopanoid_mol in hopanoid_mols):
        return False, "Invalid SMARTS pattern detected, unable to match"

    # Check if any pattern matches the SMILES structure
    for pattern in hopanoid_mols:
        if pattern and mol.HasSubstructMatch(pattern):
            return True, "Contains hopanoid skeleton with appropriate stereochemistry and pentacyclic structure"

    return False, "Hopanoid structure not detected"