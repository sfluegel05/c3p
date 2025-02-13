"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: CHEBI:35610 azole
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_azole(smiles: str):
    """
    Determines if a molecule contains an azole ring based on its SMILES string.
    An azole is a monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    The ring can also contain one or more other non-carbon atoms, such as nitrogen, sulfur, or oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an azole ring, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible azole ring pattern
    azole_pattern = Chem.MolFromSmarts("[n,s,o]1~[c,n,s,o]~[c,n,s,o]~[c,n,s,o]~[c,n,s,o]1")
    matches = mol.GetSubstructMatches(azole_pattern)
    
    if not matches:
        return False, "No five-membered ring with heteroatoms found"

    # Check each potential azole ring
    ring_info = mol.GetRingInfo()
    for match in matches:
        # Get the atoms in this ring
        ring_atoms = set()
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in match and bond.GetEndAtomIdx() in match:
                ring_atoms.add(bond.GetBeginAtomIdx())
                ring_atoms.add(bond.GetEndAtomIdx())
        
        # Check if it's a 5-membered ring
        if len(ring_atoms) != 5:
            continue
            
        # Check if it's aromatic
        is_aromatic = all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring_atoms)
        if not is_aromatic:
            continue
            
        # Check for at least one nitrogen
        has_nitrogen = any(mol.GetAtomWithIdx(atom).GetAtomicNum() == 7 for atom in ring_atoms)
        if not has_nitrogen:
            continue
            
        # Check if this ring is monocyclic (no fused rings)
        is_monocyclic = True
        for other_ring in ring_info.AtomRings():
            if set(other_ring) != ring_atoms and len(set(other_ring).intersection(ring_atoms)) > 1:
                is_monocyclic = False
                break
                
        if is_monocyclic:
            # Count heteroatoms
            heteroatoms = sum(1 for atom in ring_atoms if mol.GetAtomWithIdx(atom).GetAtomicNum() not in [6,1])
            if heteroatoms > 1:
                return True, "Contains a five-membered aromatic ring with nitrogen and other heteroatoms"
            else:
                return True, "Contains a five-membered aromatic ring with nitrogen"

    return False, "No valid azole ring found"