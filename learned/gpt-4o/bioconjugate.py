"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least two biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide bond patterns (typically featuring amide bonds)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C@H](C)C(=O)")  # Simplification: looks for alpha-amino acid amidation
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return True, "Contains peptide bond(s), indicative of protein/peptide bioconjugate"

    # Check for disulfide bond patterns
    disulfide_pattern = Chem.MolFromSmarts("S-S")
    if mol.HasSubstructMatch(disulfide_pattern):
        return True, "Contains disulfide bond(s), indicative of bioconjugate involving cysteine/dithiol"

    # Check for nucleotide-like structures linked together (common in coenzyme A type bioconjugates)
    cofactor_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C@[O,P]")  # Simplified nucleotide linkage
    if mol.HasSubstructMatch(cofactor_pattern):
        return True, "Contains cofactor-like nucleotide structure"

    # More complex logic might be needed to cover diverse examples encountered in real applications
    # Return None, None for inconclusive results based on available patterns
    return (None, None)