"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    Germacranolides are a subclass of sesquiterpene lactones typically having
    a multi-ring system including a cyclodecadiene and a lactone group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check cyclic ester (lactone); 5 or 6-membered rings are typical in germacranolides
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1OC=CC=CC=CC1"),  # generic lactone within a polycyclic framework
        Chem.MolFromSmarts("O=C1COC=CC=CC1")  # variation for different stereochemistry or substitution
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns):
        return False, "Lactone group not found"
    
    # Ensure some form of a decane or fused multi-ring system is present
    ring_info = mol.GetRingInfo()
    has_large_ring = any(len(ring) >= 10 for ring in ring_info.AtomRings())
    if not has_large_ring:
        return False, "No large ring system characteristic of germacranolides found"

    # Check multiple different double bonds present in the structure, indicating poly-unsaturation
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if num_double_bonds < 2:
        return False, "Insufficient double bonds, typical germacranolides are polyunsaturated"
    
    # Passed all checks: likely a germacranolide
    return True, "Matches typical structure of a germacranolide"