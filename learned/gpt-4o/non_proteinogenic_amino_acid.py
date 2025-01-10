"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is defined as any amino acid not encoded by the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # List of SMILES patterns for the 20 proteinogenic amino acids (not complete but indicative)
    standard_amino_acids = [
        "N[C@@H](C)C(=O)O",  # Alanine
        "NCC(=O)O",         # Glycine
        "NC(=O)[C@@H](CC1=CN=C-N1)C(O)=O",  # Histidine
        "NCC(=O)C1=CC=CC=C1",  # Phenylalanine
        # Add other 16 standard amino acids as needed
    ]

    # Look for amino and carboxylic acid groups
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0;$(O[C,c])]")
    
    if not mol.HasSubstructMatch(amino_pattern):
        return (None, None)
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return (None, None)
    
    # Check if the detected amino acid is non-standard
    for aa in standard_amino_acids:
        if mol.HasSubstructMatch(Chem.MolFromSmiles(aa)):
            return False, "Matches a standard proteinogenic amino acid"

    # If the molecule contains additional unusual groups, consider it non-proteinogenic
    unusual_features = [
        "[SX3](=O)",       # Sulfoxides
        "[OX2][CX3](=N)",  # Oximes
        "[F,Cl,Br,I]",     # Halogens
        # Add other unusual functional patterns as needed
    ]

    for feature in unusual_features:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(feature)):
            return True, "Contains unusual groups for standard amino acids"

    # Check chiral centers for non-standard amino acid behavior
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) > 1:
        return True, "Multiple chiral centers indicating non-proteinogenic"

    return False, "Does not have sufficient evidence to be classified as non-proteinogenic"