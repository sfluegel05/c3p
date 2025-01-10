"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of amino acid-like structures (common in peptides)
    amino_acid_pattern = Chem.MolFromSmarts("C([C@@H](N)C(=O)O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    if len(amino_acid_matches) < 2:
        return False, "Less than two amino acid-like structures found"
    
    # Check for certain linking groups commonly found in bioconjugates, e.g., disulfide bonds, thioethers
    disulfide_pattern = Chem.MolFromSmarts("SS")
    thioether_pattern = Chem.MolFromSmarts("CSC")
    
    disulfide_matches = mol.HasSubstructMatch(disulfide_pattern)
    thioether_matches = mol.HasSubstructMatch(thioether_pattern)
    
    if not disulfide_matches and not thioether_matches:
        return False, "No common linking groups found (disulfide or thioether)"

    # Count rotatable bonds to verify the potential linkages and flexibility
    n_rotatable = AllChem.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds to suggest linked bioconjugate"

    return True, "Contains multiple linked amino acid-like structures, indicative of a bioconjugate"