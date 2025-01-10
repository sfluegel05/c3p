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

    # Generalized Isoalloxazine core SMARTS pattern
    # The pattern represents the fused ring system in flavins, including potential variations.
    isoalloxazine_pattern = Chem.MolFromSmarts("c1cc2nc3c(nc(=O)n(c4ccc(cc4)c3=O)c2o1)")

    # Check for isoalloxazine pattern in molecule
    if mol.HasSubstructMatch(isoalloxazine_pattern):
        return True, "Contains isoalloxazine ring system characteristic of flavins"
    else:
        return False, "Does not contain isoalloxazine ring system"