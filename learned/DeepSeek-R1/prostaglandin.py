"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:15554 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are derived from prostanoic acid, characterized by a cyclopentane ring
    with two side chains, one ending in a carboxylic acid group (or derivative).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for presence of a five-membered ring
    ring_info = mol.GetRingInfo()
    five_membered_rings = [ring for ring in ring_info.AtomRings() if len(ring) == 5]
    if not five_membered_rings:
        return False, "No five-membered ring found"
    
    # Check for carboxylic acid or ester group
    acid_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    if not mol.HasSubstructMatch(acid_ester_pattern):
        return False, "No carboxylic acid/ester group"
    
    # Verify the acid/ester is attached to the five-membered ring
    acid_ester_matches = mol.GetSubstructMatches(acid_ester_pattern)
    ring_atoms = set(five_membered_rings[0])  # Consider first five-membered ring
    attached = False
    for match in acid_ester_matches:
        acid_carbon = match[0]
        # Check if acid carbon is connected to the ring via a chain
        for atom in Chem.GetShortestPath(mol, acid_carbon, ring_atoms.pop()):
            if atom in ring_atoms:
                attached = True
                break
        if attached:
            break
    if not attached:
        return False, "Carboxylic acid/ester not attached to the ring"
    
    # Check for at least two substituent chains on the ring
    substituents = 0
    for atom_idx in five_membered_rings[0]:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in five_membered_rings[0]:
                substituents += 1
                break  # Count each atom once per substituent
    if substituents < 2:
        return False, "Insufficient substituent chains on the ring"
    
    # Basic check for prostanoic acid-like structure (at least 20 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Only {carbon_count} carbons, typical prostaglandins have ≥20"
    
    return True, "Five-membered ring with carboxylic acid/ester and substituent chains"