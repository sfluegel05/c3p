"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate typically features peptide-like structures bonded
    to other types of chemical entities, often with sulfur or nitrogen elements involved.

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
    
    # SMARTS pattern to identify a broad range of peptide-like structures
    # including Cysteine & Glutathione based structures
    peptide_like_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)O)C")
    
    # SMARTS pattern for sulfur linkage, often seen in conjugation
    sulfur_conjugate_pattern = Chem.MolFromSmarts("CSC")
    
    # Check for peptide-like structures 
    if mol.HasSubstructMatch(peptide_like_pattern):
        # Also look for sulfur linkages indicative of bioconjugation
        if mol.HasSubstructMatch(sulfur_conjugate_pattern):
            return True, "Contains peptide-like structure with sulfur-linked conjugation point"
        else:
            # Additionally, look for extensive nitrogen containing structures.
            extended_nitrogen_patterns = Chem.MolFromSmarts("N=CN")
            if mol.HasSubstructMatch(extended_nitrogen_patterns):
                return True, "Contains peptide-like structure with extensive nitrogen-based conjugation chemistry"

    # If neither pattern was satisfied
    return False, "Does not match known bioconjugate structural features"