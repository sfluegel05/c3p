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

    # Sphingosine backbone pattern (typically d18:1 represents 18 carbons with one double bond)
    sphingodiene = Chem.MolFromSmarts("[C@H](O)[C@H](NC(=O)[C])[C](C=C)[C]")  # Basic pattern for sphingodiene backbone
    if not mol.HasSubstructMatch(sphingodiene):
        return False, "No sphingodiene backbone found"
    
    # Phosphocholine group
    phosphocholine_group = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")  # More inclusive pattern
    if not mol.HasSubstructMatch(phosphocholine_group):
        return False, "No phosphocholine group found"

    # Amide linkage with flexible fatty acid chain
    amide_linkage = Chem.MolFromSmarts("[NX3][C](=O)[C]")  # More general amide bond pattern
    if not mol.HasSubstructMatch(amide_linkage):
        return False, "No amide linkage found"

    return True, "Matches key patterns for sphingomyelin d18:1"