"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for carboxylic acid and amino group
    carboxylic_acid_smarts = "[CX3](=O)[OX2H1]"  # Carboxylic acid group
    amino_group_smarts = "[NX3;H2,H1;!$(N=*);!$([N+])]"  # Primary or secondary amine
    
    # Create pattern molecules
    carboxylic_acid_pattern = Chem.MolFromSmarts(carboxylic_acid_smarts)
    amino_group_pattern = Chem.MolFromSmarts(amino_group_smarts)
    
    # Search for carboxylic acid group
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    else:
        num_carboxylic_acids = len(carboxylic_acid_matches)
    
    # Search for amino group
    amino_group_matches = mol.GetSubstructMatches(amino_group_pattern)
    if not amino_group_matches:
        return False, "No amino group found"
    else:
        num_amino_groups = len(amino_group_matches)
    
    return True, f"Contains {num_carboxylic_acids} carboxylic acid group(s) and {num_amino_groups} amino group(s)"