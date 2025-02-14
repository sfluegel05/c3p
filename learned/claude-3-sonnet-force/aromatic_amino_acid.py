"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: CHEBI:33567 aromatic amino acid
An amino acid whose structure includes an aromatic ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for amino acid substructure
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)[OX2H0-]")
    amino_acid_match = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_match:
        return False, "No amino acid substructure found"
    
    # Check for aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    aromatic_ring_match = mol.GetSubstructMatches(aromatic_ring_pattern)
    if not aromatic_ring_match:
        return False, "No aromatic ring found"
    
    return True, "Contains both an amino acid substructure and an aromatic ring"