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

    # Refined SMARTS for core hopanoid skeleton, including stereochemistry
    hopanoid_smarts = "[C@]12CC[C@@H]3C[C@H](C[C@@H]4[C@]([C@H](C[C@@H]5[C@]4(CC[C@@]5(C3)C)C)(C)C)C2)C1"
    
    # Convert SMARTS into RDKit Mol object
    try:
        hopanoid_pattern = Chem.MolFromSmarts(hopanoid_smarts)
    except Exception as e:
        return None, f"Error in SMARTS pattern conversion: {e}"
    
    if hopanoid_pattern is None:
        return False, "Invalid SMARTS pattern detected, unable to match"

    # Check if pattern matches the SMILES structure
    if mol.HasSubstructMatch(hopanoid_pattern):
        return True, "Contains hopanoid skeleton with appropriate stereochemistry and pentacyclic structure"

    # Further verifications could be added here
    return False, "Hopanoid structure not detected"