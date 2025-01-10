"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    Flavins are characterized by the isoalloxazine ring system, which
    consists of a tricyclic ring system with two keto groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised Isoalloxazine pattern for flavins
    # This broader pattern should capture common structural motifs in the isoalloxazine system.
    isoalloxazine_pattern = Chem.MolFromSmarts("c1cc2nc3c(nc(=O)[nH]c3=O)c2nc1")

    # Check for isoalloxazine pattern in the molecule
    if mol.HasSubstructMatch(isoalloxazine_pattern):
        return True, "Contains isoalloxazine ring system characteristic of flavins"
    else:
        return False, "Does not contain isoalloxazine ring system"