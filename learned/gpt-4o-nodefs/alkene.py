"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene classically contains only carbon and hydrogen atoms, and at least one carbon-carbon double bond in an open-chain structure without highly reactive or complex substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carbon-carbon double bond pattern (C=C)
    alkene_pattern = Chem.MolFromSmarts("C=C")
    
    if not mol.HasSubstructMatch(alkene_pattern):
        return False, "No C=C double bond found"

    # Check for presence of heteroatoms or complex functional groups
    heteroatom_pattern = Chem.MolFromSmarts("[!C&H]")
    if mol.HasSubstructMatch(heteroatom_pattern):
        return False, "Contains heteroatom or non-carbon/hydrogen groups"

    # Verify not a small ring structure
    ring_info = mol.GetRingInfo()
    if any(len(ring) < 8 for ring in ring_info.AtomRings()):
        return False, "Contains a ring with fewer than 8 carbons"

    return True, "Contains at least one C=C double bond in suitable alkene structure"