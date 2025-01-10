"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is defined as three amino-acid residues connected by peptide linkages (amide bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for at least two amide bond patterns (N-C(=O))
    amide_bond_pattern = Chem.MolFromSmarts("N[C;R0]=O")
    amide_matches = mol.GetSubstructMatches(amide_bond_pattern)

    if len(amide_matches) < 2:
        return False, f"Less than 2 amide bonds found"

    # Check for segments connected by these amides - implies 3 residues
    # A tripeptide should only have at least 2 peptide bonds/amide-linkages i.e not including side reaction amide formations.
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)N")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)

    # Check if we effectively have two N-C(=O) spots indicative of linkages and one terminal amine or carboxyl group.
    if len(carbonyl_matches) != 2 and not (len(amide_matches) > 2):
        return False, f"Found {len(carbonyl_matches)} peptide-linkages, indicating different structure"

    # Check for the continuity of the peptide bond by assessing the sequence of N-(C=O)-C across spans
    # (simple linear peptide backbone consists of consecutive amide bonds in linear fashion)
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[NX3]")  # pattern connecting the backbone structure
    peptides = mol.GetSubstructMatches(peptide_pattern)
    if len(peptides) < 2:
        return False, "Does not have a backbone structure characteristic of peptides"

    return True, "Contains three amino-acid residues connected by peptide linkages"