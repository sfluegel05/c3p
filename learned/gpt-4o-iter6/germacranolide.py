"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a germacranolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised and expanded SMARTS patterns for germacranolide
    # Focus on 10-membered ring with lactone and sesquiterpene structural features
    germacranolide_patterns = [
        Chem.MolFromSmarts("O1[C@]2([C@](C3=CCC=C4C[C@@H](O)CC4=C3)(C(C1=O)=C)C)[H]"), # Lactone in 10-membered ring
        Chem.MolFromSmarts("C=CC1=C(C)CCC2CC(O2)C1=O"), # Germacrane core pattern
        Chem.MolFromSmarts("C1CCC2C(CC2C1)C3=CC=CC(=O)O3"), # Macrocyclic core with lactone
        # Other patterns focusing on lactones and common functional groups in sesquiterpenes
        Chem.MolFromSmarts("O=C1OC2CCC(C3=C(C=CC(=C3)C)C4=CC=CC(=O)C4)C1(C2)"), # Double-lactone & cyclic
    ]

    for pattern in germacranolide_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains germacranolide-like structure matching known patterns"

    return False, "No identifiable germacranolide-like structure using refined patterns"