"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is defined as three amino-acid residues connected by two peptide linkages (amide bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create amide bond pattern: N-C(=O)-C-N
    peptide_backbone_pattern = Chem.MolFromSmarts("N[*]-C(=O)-N[*]-C(=O)-N[*]")
    if mol.HasSubstructMatch(peptide_backbone_pattern):
        return True, "Contains three amino-acid residues connected by peptide linkages"
    
    # If no matches found, it is invalid
    return False, "Does not match peptide backbone structure"