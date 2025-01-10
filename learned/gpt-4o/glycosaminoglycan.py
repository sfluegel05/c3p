"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is characterized by polysaccharides containing substantial proportions of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more specific sugar and aminosugar patterns for glycosaminoglycans
    # Glycosaminoglycans usually consist of repeating disaccharide units containing amino sugars
    # Example of disaccharide repeat unit could be ([C@@H]1[C@@H](C(O)C(O)[C@H](O)[C@H](CO)O1)N)

    # Pattern for an aminosugar (e.g., N-acetylglucosamine)
    aminosugar_pattern = Chem.MolFromSmarts("[C@@H]([C@H](O)[C@H](O)[C@H](CO)O)N")
    if not mol.HasSubstructMatch(aminosugar_pattern):
        return False, "No aminomonosaccharide residues found"

    # Check for multiple occurrence implying a polysaccharide
    # Assuming at least 3 repeats to classify as substantial
    aminomonosaccharide_matches = mol.GetSubstructMatches(aminosugar_pattern)
    if len(aminomonosaccharide_matches) < 3:
        return False, "Does not meet the threshold for aminomonosaccharide residues"

    return True, "Contains polysaccharide chain with substantial aminomonosaccharide residues"