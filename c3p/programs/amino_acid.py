"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as having at least one carboxylic acid group
    and at least one amino group within the same molecule, excluding polymeric chains.

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
    
    # Pattern for an amino group
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")
    
    # Pattern for a carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    
    # Check for both amino and carboxyl patterns within the molecule
    has_amino = mol.HasSubstructMatch(amino_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    
    if has_amino and has_carboxyl:
        # Ensure it is not part of a known peptide or polymeric structure by looking for peptide linkages
        peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N[C]")
        if not mol.HasSubstructMatch(peptide_bond_pattern):
            return True, "Contains both amino and carboxylic acid groups without peptide linkages"
        else:
            return False, "Contains peptide bond indicative of polymers or peptide sequences"
    
    return False, "Missing amino or carboxylic acid group"