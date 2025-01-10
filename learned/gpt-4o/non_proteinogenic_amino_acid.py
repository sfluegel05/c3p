"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is defined as any amino acid not encoded by the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Comprehensive list of SMILES for proteinogenic amino acids (using SMARTS for flexibility)
    standard_amino_acids_smarts = [
        "[NX3;H2,H1;!$(NC=O)][CX4;H1,H0](C)[C](=O)[O;H1,H0]",  # Alanine
        "[NX3;H2,H1;!$(NC=O)]CC(=O)[O;H1,H0]",  # Glycine
        "[NX3;H2,H1;!$(NC=O)][CX4;H1,H0](CC1=CN=C-N1)[C](=O)[O;H1,H0]",  # Histidine
        "[NX3;H2,H1;!$(NC=O)][CX4;H1,H0](C)C1=CC=C(C=C1)[C](=O)[O;H1,H0]",  # Phenylalanine
        # Include more SMARTS for other standard amino acids
    ]

    # Look for amino and carboxylic acid groups
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0;$(O[C,c])]")

    if not mol.HasSubstructMatch(amino_pattern):
        return False, "Missing amino group"
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Missing carboxylic acid group"
    
    # Check for non-standard features
    unusual_features_smarts = [
        "[SX3](=O)",       # Sulfoxides
        "[OX2][CX3](=N)",  # Oximes
        "[F,Cl,Br,I]",     # Halogens
        "[SX2]C",          # Thioethers (sulfur bonded to carbon)
        "[Se]"             # Selenium presence
        # Add more unusual features typical in non-standard amino acids
    ]

    for feature in unusual_features_smarts:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(feature)):
            return True, "Contains unusual groups indicating non-standard"

    # If the molecule matches any of the standard patterns, return false
    for aa_pattern in standard_amino_acids_smarts:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(aa_pattern)):
            return False, "Matches a standard proteinogenic amino acid"

    # Check for additional structural complexity
    if Chem.FindMolChiralCenters(mol, includeUnassigned=True):
        return True, "Contains multiple chiral centers or unusual complexity"

    return False, "Does not match criteria for non-proteinogenic amino acid"