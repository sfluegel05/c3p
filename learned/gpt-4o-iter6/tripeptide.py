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
    
    # Define the pattern for alpha-carbon backbone N-C-C next to amide bonds
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)N")

    # Find all alpha-carbon backbone matches
    backbone_matches = mol.GetSubstructMatches(alpha_carbon_pattern)

    if len(backbone_matches) != 2:
        return (False, f"Backbone does not correspond to three contiguous residues; found: {len(backbone_matches)}")
    
    # Count peptide bonds (amide bonds)
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[NX3]")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Tripeptide must have exactly two peptide linkages characteristic of its backbone
    if len(peptide_bonds) != 2:
        return (False, f"Does not have the exact backbone characteristic of a tripeptide: {len(peptide_bonds)} found")
    
    # Check for termini or cyclic structure
    amine_pattern = Chem.MolFromSmarts("[NH2]")
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    terminal_amine = mol.HasSubstructMatch(amine_pattern)
    terminal_carboxyl = mol.HasSubstructMatch(carboxylic_pattern)
    
    # Accept cyclic or capped structures that still have 3 residues and 2 peptide bonds
    if not (terminal_amine or terminal_carboxyl or len(peptide_bonds) == 2):
        return (False, "No valid peptide termini found.")
    
    return (True, "Contains three amino-acid residues connected by peptide linkages")