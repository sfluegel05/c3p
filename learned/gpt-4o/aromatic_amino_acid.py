"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid should have at least one aromatic ring and a primary amino acid structure, 
    characterized by an amino group, a central carbon, and a carboxylic group, without additional peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for detecting any aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("a1aaaaa1")  # Represents a generic six-membered aromatic ring
    has_aromatic_ring = mol.HasSubstructMatch(aromatic_ring_pattern)
    
    # Define a SMARTS pattern for a basic amino acid backbone structure
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H2][CX4H]([CX4])([CX3](=O)[OX1H])")  # Amino group - alpha carbon - carboxylic group
    has_amino_acid_backbone = mol.HasSubstructMatch(amino_acid_pattern)
    
    # Ensure the amino acid backbone appears only once to rule out peptide chains
    is_single_amino_acid = len(mol.GetSubstructMatches(amino_acid_pattern)) == 1

    # Classification logic
    if has_aromatic_ring and has_amino_acid_backbone and is_single_amino_acid:
        return True, "Molecule is an aromatic amino acid with an aromatic ring and a simple amino acid backbone"
    elif not has_aromatic_ring:
        return False, "Missing aromatic ring"
    elif not has_amino_acid_backbone:
        return False, "Missing amino acid structure or contains peptide bonds"

    return False, "Does not match criteria for an aromatic amino acid"