"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton, featuring a pentacyclic structure with specific stereochemistry.

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
        '[C@@H]1C[C@@H]2[C@@]3(CC[C@@H]4[C@]5(C[C@@H](CC5)C[C@@H]4C[C@]3(C)CC2)C)C1',
        '[C@H]1[C@@H]2CC[C@@]3(C)[C@]([C@@H]4C[C@@H](CC4)C[C@@H]3C)[C@H]2CC1',
        '[C@@H]1CC[C@@]2([C@H]([C@@H]3C[C@@H](CC3)[C@@]2(C)CC1)C)C',
        # More hopanoid patterns can be added here
    ]

    # Convert SMARTS to RDKit Mol objects
    try:
        hopanoid_mols = [Chem.MolFromSmarts(pattern) for pattern in hopanoid_patterns]
    except Exception as e:
        return None, f"Error in SMARTS pattern conversion: {e}"
    
    if any(hopanoid_mol is None for hopanoid_mol in hopanoid_mols):
        return False, "Invalid SMARTS pattern detected, unable to match"

    # Check for matches with any of the defined patterns
    for pattern in hopanoid_mols:
        if pattern and mol.HasSubstructMatch(pattern):
            return True, "Contains hopanoid skeleton with appropriate stereochemistry and pentacyclic structure"

    return False, "Hopanoid structure not detected"