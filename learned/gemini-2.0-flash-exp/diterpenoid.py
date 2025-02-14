"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids have a C20 skeleton (or modified) and are built from isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for ring systems
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    
    if n_rings < 2:
        return False, f"Diterpenoids should have at least 2 rings, this one has {n_rings}"

    # Check for presence of multiple C5 units. Count C atoms with 1 or 2 H, or methyls.
    # [CX4H2,CX4H1,CX4H3]
    c5_pattern = Chem.MolFromSmarts("[CX4H2,CX4H1,CX4H3]~[CX4H2,CX4H1,CX4H3]~[CX4H2,CX4H1,CX4H3]~[CX4H2,CX4H1,CX4H3]~[CX4H2,CX4H1,CX4H3]")
    c5_matches = mol.GetSubstructMatches(c5_pattern)

    if len(c5_matches) < 1 : #Requires at least one C5 unit
         return False, "Too few connected isoprene units detected (less than 1)"

    # Complex ring check to avoid simple bicyclic systems, require a fused ring system.
    # Look for rings where two rings share 2 or more atoms.
    fused_ring_pattern = Chem.MolFromSmarts("*1~*~*~*~*~*1*2~*~*~*~*~*2")
    fused_ring_matches = mol.GetSubstructMatches(fused_ring_pattern)

    if len(fused_ring_matches) < 1:
        return False, "Too few fused rings"


    return True, "Likely a diterpenoid based on ring system and multiple isoprene units"