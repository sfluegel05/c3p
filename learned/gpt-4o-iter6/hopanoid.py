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

    # Hopanoid hopane skeleton pattern - five fused rings with specific stereochemistry
    hopanoid_pattern = Chem.MolFromSmarts('[C@]12CC[C@@H]3[C@]4(C)CC[C@@]5([C@H]4CC[C@]35C)C2C1')
    if not mol.HasSubstructMatch(hopanoid_pattern):
        return False, "Hopanoid hopane skeleton not found"

    return True, "Contains hopanoid hopane skeleton based on pentacyclic triterpenoid structure"