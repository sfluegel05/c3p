"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide comprises a peptide moiety and an attached lipid (long hydrocarbon chain).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify peptide bonds - expand amide bond pattern for flexibility
    peptide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(peptide_pattern):
        return False, "No peptide bonds found"
    
    # Flexible identification for long hydrocarbon chains
    # Using pattern for long continuous aliphatic chains (at least 8 carbon atoms)
    lipid_pattern = Chem.MolFromSmarts("CCCCCCCC")  # Simple form of a long chain pattern
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No long hydrocarbon chains found"

    return True, "Contains both peptide bonds and long hydrocarbon chains, characteristic of lipopeptides"