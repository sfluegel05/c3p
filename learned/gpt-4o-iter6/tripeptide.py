"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of three amino-acid residues connected by peptide linkages (amide bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Define the pattern for amide bonds denoting peptide linkage
    amide_bond_pattern = Chem.MolFromSmarts("N[C](=O)")
    
    # Find all amide bond matches
    amide_matches = mol.GetSubstructMatches(amide_bond_pattern)

    if len(amide_matches) < 2:
        return (False, f"Less than 2 amide bonds found; found: {len(amide_matches)}")
    
    # Count sequence of N-C(=O)-C (should be two for precise tripeptide with three residues)
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[NX3]")
    peptides = mol.GetSubstructMatches(peptide_pattern)
    
    # Tripeptide must have exactly two linkages characteristic of its backbone
    if len(peptides) != 2:
        return (False, f"Does not have the exact backbone characteristic of a tripeptide: {len(peptides)} found")
    
    # Check for expected termini (either NH2 and COOH or other cyclic ends)
    amine_pattern = Chem.MolFromSmarts("[NH2]")
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    terminal_amine = mol.HasSubstructMatch(amine_pattern)
    terminal_carboxyl = mol.HasSubstructMatch(carboxylic_pattern)
    
    if not (terminal_amine or terminal_carboxyl):
        return (False, "No valid peptide termini found")
    
    return (True, "Contains three amino-acid residues connected by peptide linkages")