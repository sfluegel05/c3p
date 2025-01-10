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
        bool: True if the molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # List of SMARTS for all standard amino acids (for illustration purposes, cover more in practice)
    standard_amino_acids_smarts = [
        "[NX3;H2,H1;!$(NC=O)][CX4;H1,H0](C(O)=O)",  # Alanine example
        # Add SMARTS patterns for all 20 standard amino acids
    ]

    # Check for presence of amino and carboxylic acid groups
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0;!$(NC)]")

    if not mol.HasSubstructMatch(amino_pattern):
        return False, "Missing amino group"
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Missing carboxylic acid group"
    
    # Detect unusual features typical in non-standard amino acids
    unusual_features_smarts = [
        "[SX3](=O)",       # Sulfoxides
        "[OX2][CX3](=N)",  # Oximes
        "[F,Cl,Br,I]",     # Halogens
        "[SX2]C",          # Thioethers
        "[Se]",            # Selenium
        "[PH]O"            # Phosphorus with oxygen
        # Extend with more patterns
    ]

    for feature in unusual_features_smarts:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(feature)):
            return True, "Contains unusual groups indicating non-standard"

    # Check for extended complexity beyond proteinogenic structures (rings, branches)
    complexity_patterns = [
        "[C](=[O,N])O[C,H]",  # Ester linkages uncommon in standard
        "[N][C](=C)"          # An amine connected to a vinyl group
        # Add more complexity variations
    ]

    for pattern in complexity_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, "Contains structural complexity indicating non-standard"

    for aa_pattern in standard_amino_acids_smarts:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(aa_pattern)):
            return False, "Matches a standard proteinogenic amino acid"

    return False, "Does not match criteria for non-proteinogenic amino acid"