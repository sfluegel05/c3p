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

    # Define a SMARTS pattern for generic aminomonosaccharides (looking for -NH or -NH2 on sugars)
    aminomonosaccharide_patterns = [
        Chem.MolFromSmarts("[C@H]1([NH2])[C@@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])[O@H]1"),  # Example of glucosamine
        Chem.MolFromSmarts("[C@H]1([NH][C@H](C=O)*)[C@@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])[O@H]1")  # Example of N-acetylglucosamine
    ]
    
    total_aminomonosaccharide_matches = 0
    for pattern in aminomonosaccharide_patterns:
        matches = mol.GetSubstructMatches(pattern)
        total_aminomonosaccharide_matches += len(matches)

    # Estimate the total number of sugar-like units by counting oxygens (rough heuristic)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if total_aminomonosaccharide_matches > 0 and total_aminomonosaccharide_matches / (o_count / 2) > 0.2:
        # Assuming a glycosaminoglycan would have aminomonosaccharide in a substantial proportion of units
        return True, f"Glycosaminoglycan detected, substantial aminomonosaccharide residues found: count is {total_aminomonosaccharide_matches}"
    
    return False, "No substantial proportion of aminomonosaccharide residues found"