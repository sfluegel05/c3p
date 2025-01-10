"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans are polysaccharides with a significant proportion of aminomonosaccharide residues.

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

    # Define patterns for possible aminomonosaccharides, including variants
    aminomonosaccharide_patterns = [
        Chem.MolFromSmarts("[C@H]1([NH2])[C@@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])[O@H]1"),  # Glucosamine type
        Chem.MolFromSmarts("[C@H]1([NH2])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]([OH])[O@H]1"),  # Galactosamine type
        Chem.MolFromSmarts("[C@H]1([NH][C@H](C=O)*)[C@@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])[O@H]1"),  # N-acetylglucosamine type
        Chem.MolFromSmarts("[C@@H]1([C@H](O)[C@H](O)[C@@H](N)[C@@H]1O)"),  # General aminomonosaccharide
    ]
    
    total_aminomonosaccharide_matches = 0
    for pattern in aminomonosaccharide_patterns:
        matches = mol.GetSubstructMatches(pattern)
        total_aminomonosaccharide_matches += len(matches)

    # Calculate a better estimation for sugar units
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Glycosaminoglycans are largely composed of sugars, let's assume ~3 oxygens per sugar for better specificity
    sugar_units_estimation = max(1, oxygen_count // 3)

    # Set threshold for classification, assuming at least 30% units should be amino sugars
    if total_aminomonosaccharide_matches / sugar_units_estimation > 0.3:
        return True, f"Glycosaminoglycan identified. Aminomonosaccharide count: {total_aminomonosaccharide_matches}"

    return False, "Insufficient proportion of aminomonosaccharide residues found"