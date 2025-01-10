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
    
    # Identify peptide bonds
    peptide_pattern = Chem.MolFromSmarts("C(=O)N")  # Amide bond pattern
    if not mol.HasSubstructMatch(peptide_pattern):
        return False, "No peptide bonds found"
    
    # Improved pattern for flexible identification of long hydrocarbon chains
    # Allow for some branching - target a minimum of 8 carbon atoms overall
    lipid_pattern = Chem.MolFromSmarts("C~C~C~C~C~C~C~C")  # Flexible pattern with connectivity
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No sufficient long hydrocarbon chains found"

    return True, "Contains both peptide bonds and sufficient long hydrocarbon chains, characteristic of lipopeptides"