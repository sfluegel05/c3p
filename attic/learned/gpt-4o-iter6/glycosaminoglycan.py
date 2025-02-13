"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is a polysaccharide containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for aminomonosaccharides
    aminomonosaccharide_pattern = Chem.MolFromSmarts("[$([NX3;H2,H1;!$(NC=O)])]-[*]-[*]-[*]-[*]")  # N attached to a basic carbohydrate substructure

    # Find all matches for the aminomonosaccharide substructure
    matches = mol.GetSubstructMatches(aminomonosaccharide_pattern)
    num_aminomonosaccharide = len(matches)

    # Use a basic heuristic to determine if the proportion is substantial
    if num_aminomonosaccharide > 0:
        return True, f"Contains aminomonosaccharide residues, count: {num_aminomonosaccharide}"
    
    return False, "No substantial proportion of aminomonosaccharide residues found"