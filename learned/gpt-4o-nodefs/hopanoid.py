"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    Hopanoids are a class of pentacyclic triterpenoids with hopane backbones.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated patterns capturing a broad set of characteristics for hopanoids
    hopanoid_patterns = [
        # Basic hopane structure with possible methyl groups
        Chem.MolFromSmarts("C1CC[C@]2(C)[C@H](CC[C@@H]3[C@@]4(CC[C@H]5[C@@]6(C)CCCC(C)(C)[C@@H]6CC[C@@H]5C)C)[C@@H](C4)CC3C2)C1"),
        # Consider variations where cyclopentane motifs might differ
        Chem.MolFromSmarts("C1CCC2(CC3CC4(C)C5CCC5C4CC3C2C1)C"),
        # Allow for functionalized hopanoids (with flexible oxygen or nitrogens)
        Chem.MolFromSmarts("C1[C@H]2C[C@]3(CC[C@@H]4[C@@]5(C(C)(C)C)CCCC(C@@H]5CC[C@@H]4C3C2)[C@@H](C1)[NX2,O]"),
    ]

    for hopanoid_pattern in hopanoid_patterns:
        if mol.HasSubstructMatch(hopanoid_pattern):
            return True, "Contains hopanoid-like backbone"

    return False, "No hopanoid backbone found"

# This code includes broader SMARTS patterns that cover hopanoid structures, considering common variations
# in functional groups, stereochemistry, and ring structures.