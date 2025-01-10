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
    
    # Check for a 5- or 6-membered lactone ring with considering diverse linkages
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1COC=C1"),  # 5-membered lactone with unsaturation
        Chem.MolFromSmarts("O=C1CCCCO1"),  # 6-membered lactone without unsaturation
        Chem.MolFromSmarts("O=C1COCC1"),   # 5-membered lactone variation
        Chem.MolFromSmarts("O1CC=COC1=O")  # 6-membered lactone with double bonds
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns):
        return False, "Lactone group not found"
    
    # Ensure some form of a decane or fused multi-ring system is present
    ring_info = mol.GetRingInfo()
    has_large_ring = any(len(ring) >= 8 for ring in ring_info.AtomRings())  # Modified threshold to be more inclusive
    if not has_large_ring:
        return False, "No large ring system characteristic of germacranolides found"

    # Check for multiple double bonds present indicating poly-unsaturation
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if num_double_bonds < 2:
        return False, "Insufficient double bonds, typical germacranolides are polyunsaturated"
    
    # Consider stereochemistry as it might be crucial in classification
    # However, it's typically complicated without visual assessment, so this step is non-mandatory unless specific needs arise

    return True, "Matches typical structure of a germacranolide"