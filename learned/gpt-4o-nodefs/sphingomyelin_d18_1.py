"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    Sphingomyelin d18:1 has a specific sphingosine backbone, a phosphocholine group, 
    and amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"

    # Sphingosine backbone with specific unsaturation (trans- or cis-)
    sphingosine_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](NC(=O)[C])(/C=C/C)")  # Adjusted to match common patterns
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone with correct stereochemistry found"
    
    # Phosphocholine head group
    phosphocholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Amide linkage to a fatty acid
    amide_pattern = Chem.MolFromSmarts("NC(=O)C")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    return True, "Matches key patterns for sphingomyelin d18:1"