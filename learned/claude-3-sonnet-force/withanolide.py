"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: CHEBI:52356 withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is a C28 steroid lactone with a modified side chain forming a lactone ring and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid skeleton (4 fused rings with specific ring sizes)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    steroid_ring_sizes = [5, 6, 6, 6]
    if set(ring_sizes) != set(steroid_ring_sizes):
        return False, "Does not contain steroid skeleton (4 fused rings with specific ring sizes)"

    # Check for lactone ring
    lactone_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "Does not contain lactone ring"

    # Check for 28 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 28:
        return False, f"Does not have 28 carbon atoms, found {c_count}"

    # Check for modified side chain (not a simple alkyl chain)
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    alkyl_chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if len(alkyl_chain_matches) > 1:
        return False, "Contains simple alkyl side chain, not modified"

    # Check for lactone ring formed by the side chain
    side_chain_lactone_pattern = Chem.MolFromSmarts("[OX2]~[CX3](=[OX1])~[CX3]~[CX3]~[CX3]")
    if not mol.HasSubstructMatch(side_chain_lactone_pattern):
        return False, "Lactone ring not formed by the side chain"

    # Check for presence of substituents on the steroid core
    substituents_pattern = Chem.MolFromSmarts("[CX3]~[OX2]")
    if not mol.HasSubstructMatch(substituents_pattern):
        return False, "No substituents on the steroid core"

    # Passed all tests, classify as withanolide
    return True, "Contains steroid skeleton with modified side chain forming lactone ring and its substituted derivatives"