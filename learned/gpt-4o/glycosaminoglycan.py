"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is characterized by polysaccharides containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more generalized pattern for aminosugar
    # Aminosugars typically have an NH group attached to a sugar (ring structure with carbon and oxygen)
    aminosugar_patterns = [
        Chem.MolFromSmarts("CC(O)C1OC(CO)C(O)C1N"),  # Pattern similar to glucosamine
        Chem.MolFromSmarts("CC(O)C1OC(CO)C(N)C1O"),  # Pattern similar to galactosamine
        Chem.MolFromSmarts("CC(O)C1OC(C)C(N)C1O"),   # Considering variations
    ]
    
    # Check if any aminosugar pattern matches the structure
    aminosugar_match_found = any(mol.HasSubstructMatch(pattern) for pattern in aminosugar_patterns)
    if not aminosugar_match_found:
        return False, "No aminomonosaccharide residues found"

    # Count matching substructures as a form of quantification
    total_matches = sum(len(mol.GetSubstructMatches(pattern)) for pattern in aminosugar_patterns)
    
    # Determine if there's a substantial polymetric chain (as heuristics, more than 5 subunits)
    if total_matches < 5:
        return False, f"Found {total_matches} aminomonosaccharide matches, need more for classification"

    return True, "Polysaccharide chain with substantial aminomonosaccharide residues found"