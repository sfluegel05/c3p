"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate typically features peptide-like structures bonded to other chemical entities,
    with linkage involving sulfur, nitrogen, or other bioconjugation relevant elements.

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
    
    # SMARTS patterns to identify a broad range of peptide-like structures
    peptide_like_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)O)")
    
    # SMARTS patterns for various possible sulfur and nitrogen linkages
    # Beyond simple "CSC", consider "CSN", "SNS", as well as any O-N linkages
    sulfur_nitrogen_linkage_patterns = [
        Chem.MolFromSmarts("CSC"), 
        Chem.MolFromSmarts("CSN"), 
        Chem.MolFromSmarts("SNS"), 
        Chem.MolFromSmarts("C-N-C"),
        Chem.MolFromSmarts("C(=O)-N"),
    ]
    
    # Check for peptide-like structures
    if mol.HasSubstructMatch(peptide_like_pattern):
        for pattern in sulfur_nitrogen_linkage_patterns:
            if mol.HasSubstructMatch(pattern):
                return True, "Contains peptide-like structure with relevant conjugation point"

    return False, "Does not match known bioconjugate structural features"