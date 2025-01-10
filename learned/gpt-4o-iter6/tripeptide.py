"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is defined as three amino-acid residues connected by two peptide linkages (amide bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide bond pattern (-N-C(=O)-)
    amide_pattern = Chem.MolFromSmarts("N-C(=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 2:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 2"

    # Count carboxylic acids and primary amines to verify amino acid residues
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = len(mol.GetSubstructMatches(acid_pattern))
    
    amine_pattern = Chem.MolFromSmarts("N")
    amine_matches = len(mol.GetSubstructMatches(amine_pattern))

    if acid_matches < 3 or amine_matches < 3:
        return False, f"Need at least 3 carboxylic acids and amines for 3 residues, got {acid_matches} acids and {amine_matches} amines"
    
    return True, "Contains three amino-acid residues connected by peptide linkages"