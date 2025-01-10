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

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for aminomonosaccharides
    aminomonosaccharide_patterns = [
        Chem.MolFromSmarts("[C@H]([NH2])[C@H]([OH])[C@H]([OH])[C@H]([OH])O"),  # General aminosugar
        Chem.MolFromSmarts("C([NH][C=O])([OH])C([OH])C([OH])O"),  # Acetylated amino sugar
        Chem.MolFromSmarts("[C@H]([NH2])C([OH])([OH])C([OH])C([OH]")  # Another variant
    ]

    # Verify patterns creation, avoid patterns that return None
    for pattern in aminomonosaccharide_patterns:
        if pattern is None:
            return (None, "Failed to create a valid SMARTS pattern for an aminomonosaccharide")

    # Count matching substructures for aminomonosaccharides
    total_aminomonosaccharide_matches = 0
    for pattern in aminomonosaccharide_patterns:
        matches = mol.GetSubstructMatches(pattern)
        total_aminomonosaccharide_matches += len(matches)

    # Estimate number of sugar units from oxygen atoms (heuristic)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    sugar_units_estimation = max(1, oxygen_count // 3)  # Rough estimate of sugars

    # Determine if a significant portion are aminomonosaccharides
    if total_aminomonosaccharide_matches / sugar_units_estimation >= 0.25:
        return True, f"Glycosaminoglycan identified. Aminomonosaccharide count: {total_aminomonosaccharide_matches}"

    return False, "Insufficient proportion of aminomonosaccharide residues found"