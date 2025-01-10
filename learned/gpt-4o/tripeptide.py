"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
        
    # Check for three peptide bonds (N-C(=O)-C)
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)-C")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 2:
        return False, "Insufficient peptide bonds, need at least 2 to form a tripeptide"
    
    # Count number of peptides
    n_peptides = len(peptide_bond_matches) + 1  # one more than peptide bonds
    if n_peptides != 3:
        return False, f"Molecule has {n_peptides} amino acids, should have 3"
    
    # Ensure amino acid backbone (Cα - N - C(=O) - C backbone pattern)
    amino_acid_backbone_pattern = Chem.MolFromSmarts("[CX4H]-[NH]-C(=O)-[CX3]")
    if not mol.HasSubstructMatch(amino_acid_backbone_pattern):
        return False, "No amino acid backbone detected"
    
    return True, "Molecule is a tripeptide with three amino acids connected by peptide linkages"