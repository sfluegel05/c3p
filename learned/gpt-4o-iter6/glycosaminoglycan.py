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

    # Define extended SMARTS patterns for aminomonosaccharides
    aminomonosaccharide_patterns = [
        Chem.MolFromSmarts("[C@H]1([NH2])[C@@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])[O@H]1"),  # Glucosamine
        Chem.MolFromSmarts("[C@H]1([NH2])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])[C@H]([OH])[O@H]1"),  # Galactosamine
        Chem.MolFromSmarts("[C@H]1([NH][C@H](C=O)*)[C@@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])[O@H]1")  # N-acetylglucosamine
    ]
    
    total_aminomonosaccharide_matches = 0
    for pattern in aminomonosaccharide_patterns:
        matches = mol.GetSubstructMatches(pattern)
        total_aminomonosaccharide_matches += len(matches)

    # Estimate the total number of sugar-like units by counting oxygen atoms
    oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    sugar_units_estimation = oxygen_atoms // 2  # heuristic: two oxygens per sugar unit

    # Adjusting the percentage threshold to define substantial as 30%
    if total_aminomonosaccharide_matches > 0 and total_aminomonosaccharide_matches / sugar_units_estimation > 0.3:
        return True, f"Glycosaminoglycan identified. Aminomonosaccharide count: {total_aminomonosaccharide_matches}"

    return False, "Insufficient proportion of aminomonosaccharide residues found"