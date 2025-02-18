"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as having a carboxylic acid group and an amino group
    both connected to a potential common carbon center (often the alpha carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classifiable as an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern allowing for variation in the connection of amino and carboxyl groups
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H2,NX3H1][CX4]([CX3](=O)[OX2H1])")
    
    # This pattern attempts to recognize an amino group possibly attached to a carbon that also connects to a carboxyl group
    if mol.HasSubstructMatch(amino_acid_pattern):
        return True, "Contains amino acid backbone structure"

    # Additional check: compound must not be a peptide or other polymer but single amino acid or modified single amino acid
    # This can be refined and expanded with more detailed structural checks or multiple patterns
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Contains peptide bond indicative of polymers or peptide sequences"

    return False, "No amino acid backbone structure found"