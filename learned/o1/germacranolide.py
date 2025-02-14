"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton,
    which is a 10-membered ring system, typically with a gamma-lactone ring fused to it.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Find 10-membered rings
    ten_membered_rings = [set(ring) for ring in atom_rings if len(ring) == 10]
    if not ten_membered_rings:
        return False, "No 10-membered ring (germacrane skeleton) found"
    
    # Find 5-membered rings (potential lactone rings)
    five_membered_rings = [set(ring) for ring in atom_rings if len(ring) == 5]
    if not five_membered_rings:
        return False, "No 5-membered ring (potential lactone ring) found"
    
    # Identify lactone rings (5-membered cyclic esters)
    lactone_smarts = "[O]=[C]-O-[C;R]"
    lactone_mol = Chem.MolFromSmarts(lactone_smarts)
    lactone_matches = mol.GetSubstructMatches(lactone_mol)
    lactone_atoms = [set(match) for match in lactone_matches]
    
    lactone_rings = []
    for ring in five_membered_rings:
        for lactone in lactone_atoms:
            # Check if the lactone substructure is within the ring
            if lactone.issubset(ring):
                lactone_rings.append(ring)
                break
    if not lactone_rings:
        return False, "No lactone ring found among 5-membered rings"
    
    # Check for fused rings (shared atoms between 10-membered ring and lactone ring)
    fused = False
    for ten_ring in ten_membered_rings:
        for lactone_ring in lactone_rings:
            shared_atoms = ten_ring & lactone_ring
            if shared_atoms:
                fused = True
                break
        if fused:
            break
    if not fused:
        return False, "No fused lactone ring found with 10-membered ring"
    
    # Optionally, check for sesquiterpene backbone (15 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Too few carbons ({c_count}) for a germacranolide (sesquiterpene)"
    
    return True, "Molecule is a germacranolide (10-membered ring fused with lactone ring)"