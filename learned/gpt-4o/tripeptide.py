"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of three amino acid residues connected by peptide linkages.

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
    
    # Peptide bond pattern - typical connection between amino acids [N-C(=O)-]
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;!H0]=C")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Ensure exactly 2 peptide bonds since tripeptides consist of three residues connected by two peptide bonds
    if len(peptide_bond_matches) != 2:
        return False, f"Expected 2 peptide bonds, found {len(peptide_bond_matches)}"
    
    # Amino acid identification based on the central structure includes typical variations.
    amino_acid_varied_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[CX4]")
    backbone_matches = mol.GetSubstructMatches(amino_acid_varied_pattern)
    n_amino_acids = len(set(match[0] for match in backbone_matches))  # Unique amino group presence
    
    # Validate exactly 3 amino acids
    if n_amino_acids != 3:
        return False, f"Molecule spans {n_amino_acids} amino acids, should have 3"
    
    return True, "Molecule is a tripeptide with three amino acids connected by peptide linkages"