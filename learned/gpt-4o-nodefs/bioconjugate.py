"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate typically includes peptide-like structures connected via
    sulfur, nitrogen, phosphorus linkages or other complex bio-relevant motifs.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to identify peptide-like structures
    peptide_like_patterns = [
        Chem.MolFromSmarts("N[C@H](C)C(=O)N"),  # Common peptide bond
        Chem.MolFromSmarts("N[C@H](C(=O))"),   # Peptide with specific stereo
        Chem.MolFromSmarts("C(=O)N[C@H]")      # Reverse peptide bond
    ]
    
    # SMARTS patterns for specific conjugation points found in bioconjugates
    conjugation_point_patterns = [
        Chem.MolFromSmarts("CSC"),  # Simple sulfur linkages
        Chem.MolFromSmarts("CSN"),  # Sulfur-nitrogen
        Chem.MolFromSmarts("SSC"),  # Disulfide linkage
        Chem.MolFromSmarts("COPO"), # Phosphorylated
        Chem.MolFromSmarts("C-P-O"),# Phosphate linkages
    ]
    
    # Check for peptide-like structures
    has_peptide = any(mol.HasSubstructMatch(pattern) for pattern in peptide_like_patterns)
    if has_peptide:
        # Check for conjugation points indicative of bioconjugates
        for pattern in conjugation_point_patterns:
            if mol.HasSubstructMatch(pattern):
                return True, "Contains peptide-like structure with relevant conjugation point"

    return False, "Does not match known bioconjugate structural features"