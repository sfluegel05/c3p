"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    Bioconjugates are characterized by the presence of peptide-like structures covalently bonded
    to other chemical entities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Common patterns for amino acids such as Cysteine or Glutathione
    amino_acid_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)O)C")

    # Check for amino acid-like structures 
    if mol.HasSubstructMatch(amino_acid_pattern):
        # Check for significant number of sulfur (S) or nitrogen (N) atoms that could indicate conjugation points
        s_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'S')
        n_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N')
        
        # Set a threshold for counting sulfur and nitrogen atoms to infer potential conjugation
        if s_count >= 2 or n_count >= 3:
            return True, "Contains amino acid-like structure with potential conjugation points"
    
    return False, "Does not match known bioconjugate structural features"