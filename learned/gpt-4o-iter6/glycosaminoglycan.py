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

    # Define patterns for aminomonosaccharides with broader generalization
    aminomonosaccharide_patterns = [
        # More generic aminosugar patterns captured
        Chem.MolFromSmarts("[C@H]([NH2])C([OH])C([OH])C([OH])C([OH])O"),  # General
        Chem.MolFromSmarts("[C@H]([NH][C@H](C=O)])C([OH])C([OH])C([OH])O")  # Acetylated
    ]

    total_aminomonosaccharide_matches = 0
    for pattern in aminomonosaccharide_patterns:
        matches = mol.GetSubstructMatches(pattern)
        total_aminomonosaccharide_matches += len(matches)

    # Estimate number of sugar units from oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    sugar_units_estimation = max(1, oxygen_count // 3)  # Rough estimate of sugars

    # Modify the threshold to determine glycosaminoglycans
    # Assuming a significant portion >= 25% of the sugar units should be aminomonosaccharide
    if total_aminomonosaccharide_matches / sugar_units_estimation >= 0.25:
        return True, f"Glycosaminoglycan identified. Aminomonosaccharide count: {total_aminomonosaccharide_matches}"

    return False, "Insufficient proportion of aminomonosaccharide residues found"