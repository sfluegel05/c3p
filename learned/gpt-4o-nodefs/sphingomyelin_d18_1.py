"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    Sphingomyelin d18:1 has a d18:1 sphingosine backbone, a phosphocholine headgroup, 
    and an N-acyl chain (amide-linked fatty acid).

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

    # Pattern for d18:1 sphingosine backbone: long chain with a hydroxyl at C3, a trans double bond typically between C4 and C5
    sphingosine_backbone = Chem.MolFromSmarts("[C@H](O)[C@H](NC(=O)C)[C](C=C)CC")
    
    if not mol.HasSubstructMatch(sphingosine_backbone):
        return False, "No sphingosine d18:1 backbone found"
    
    # Phosphocholine group
    phosphocholine_group = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C") 
    if not mol.HasSubstructMatch(phosphocholine_group):
        return False, "No phosphocholine group found"

    # Amide linkage for N-acyl chain
    amide_linkage = Chem.MolFromSmarts("NC(=O)C")  # Adapted to identify the N-acyl connection
    if not mol.HasSubstructMatch(amide_linkage):
        return False, "No amide linkage found"

    return True, "Structure matches critical patterns for sphingomyelin d18:1"