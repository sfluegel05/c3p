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
    
    # Identify peptide bonds - look for amide bonds using SMARTS: C(=O)N
    peptide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(peptide_pattern):
        return False, "No peptide bonds found"
    
    # Identify long hydrocarbon chains typical of lipids
    # Long aliphatic chain pattern: CH2 chain of length 6 or greater
    lipid_pattern = Chem.MolFromSmarts("C(CCCCC)C")
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No long hydrocarbon chains found"

    return True, "Contains both peptide bonds and long hydrocarbon chains, characteristic of lipopeptides"