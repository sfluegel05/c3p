"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: CHEBI:51922 withanolide
A withanolide is any steroid lactone that is a C28 steroid with a modified side chain forming a lactone ring and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.

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

    # Look for steroid core (4 fused rings with specific ring sizes)
    ring_info = mol.GetRingInfo()
    cores = ring_info.AtomRings()
    steroid_core = False
    for core in cores:
        ring_sizes = [ring_info.RingSize(ring) for ring in core]
        if len(core) == 4 and set(ring_sizes) == {5, 6, 6, 6}:
            steroid_core = True
            break
    if not steroid_core:
        return False, "No steroid core found"

    # Look for lactone ring
    lactone_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone ring found"

    # Check for modified side chain
    side_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[OX2][CX3](=[OX1])")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if not side_chain_matches:
        return False, "No modified side chain forming a lactone ring found"

    # Check carbon count (allow some flexibility)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 26 or c_count > 32:
        return False, f"Carbon count ({c_count}) outside expected range for withanolides"

    # Check molecular weight (typically >400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for withanolide"

    return True, "Contains steroid core with modified side chain forming a lactone ring"