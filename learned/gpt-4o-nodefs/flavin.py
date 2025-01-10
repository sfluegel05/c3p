"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    Flavins are heterocyclic compounds characterized by the isoalloxazine ring system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised Isoalloxazine core SMARTS pattern
    # This pattern seeks to more accurately match the structure typically found in flavins.
    isoalloxazine_pattern = Chem.MolFromSmarts("c1cc2nc3c(c(=O)n(c4ccccc4)c3=O)c2nc1")

    # Check for isoalloxazine pattern in molecule
    if mol.HasSubstructMatch(isoalloxazine_pattern):
        return True, "Contains isoalloxazine ring system characteristic of flavins"
    else:
        return False, "Does not contain isoalloxazine ring system"