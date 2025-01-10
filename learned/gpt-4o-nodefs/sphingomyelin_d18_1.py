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

    # More accurate pattern for d18:1 sphingosine backbone: 18 carbon atoms, a hydroxyl on C1 or C3, and one double bond
    # This SMARTS pattern still requires knowledge of the specific configuration and positions, assumes C1 hydroxyl
    sphingosine_backbone = Chem.MolFromSmarts("C[C@H](O)CCCCCCCC=C[C@H](NC(=O)C)") 
    if not mol.HasSubstructMatch(sphingosine_backbone):
        return False, "No sphingosine d18:1 backbone found"
    
    # Check for phosphocholine headgroup
    phosphocholine_group = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_group):
        return False, "No phosphocholine group found"

    # Check for amide linkage indicating N-acyl chain
    amide_linkage = Chem.MolFromSmarts("NC(=O)C")
    if not mol.HasSubstructMatch(amide_linkage):
        return False, "No amide linkage found"

    return True, "Structure matches critical patterns for sphingomyelin d18:1"