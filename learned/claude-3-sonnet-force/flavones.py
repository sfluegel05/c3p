"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: CHEBI:17794 flavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone is a flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and remove explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    Chem.RemoveHs(mol)
    
    # Kekulize to handle tautomers
    Chem.Kekulize(mol)
    
    # Look for flavone skeleton pattern
    flavone_skeleton = has_flavone_skeleton(mol)
    if not flavone_skeleton:
        return False, "Does not contain the flavone skeleton"
    
    # Count aromatic rings
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings != 3:
        return False, f"Expected 3 aromatic rings, found {aromatic_rings}"
    
    # Check for oxygens
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 4:
        return False, "Too few oxygens for a flavone"
    
    # Check for hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts("OC")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    num_hydroxy = len(hydroxy_matches)
    if num_hydroxy < 1:
        return False, "No hydroxyl groups found"
    
    # Optional: Check molecular weight and Lipinski's rules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside the typical range for flavones"
    
    if not rdMolDescriptors.CalcLipinski(mol):
        return False, "Violates Lipinski's rules for drug-likeness"
    
    return True, "Contains the flavone skeleton with expected aromatic rings, oxygens, and hydroxyl groups"

def has_flavone_skeleton(mol):
    """
    Checks if the molecule contains the flavone skeleton (2-aryl-1-benzopyran-4-one or 2-arylchromen-4-one).

    Args:
        mol (Mol): RDKit molecule object

    Returns:
        bool: True if the molecule contains the flavone skeleton, False otherwise
    """
    rings = mol.GetRingInfo().AtomRings()
    for ring1, ring2 in itertools.combinations(rings, 2):
        if len(ring1) == 6 and len(ring2) == 6 and set(ring1) & set(ring2):
            # Two fused rings, check if one is benzene and the other is pyran/pyrone
            ring1_atoms = [mol.GetAtomWithIdx(idx) for idx in ring1]
            ring2_atoms = [mol.GetAtomWithIdx(idx) for idx in ring2]
            
            if is_benzene_ring(ring1_atoms) and is_pyran_pyrone_ring(ring2_atoms):
                return True
            elif is_benzene_ring(ring2_atoms) and is_pyran_pyrone_ring(ring1_atoms):
                return True
    
    return False

def is_benzene_ring(ring_atoms):
    """
    Checks if a ring is a benzene ring.

    Args:
        ring_atoms (list): List of Atom objects in the ring

    Returns:
        bool: True if the ring is a benzene ring, False otherwise
    """
    return all(atom.GetAtomicNum() == 6 and atom.GetIsAromatic() for atom in ring_atoms)

def is_pyran_pyrone_ring(ring_atoms):
    """
    Checks if a ring is a pyran or pyrone ring with a carbonyl group at position 4.

    Args:
        ring_atoms (list): List of Atom objects in the ring

    Returns:
        bool: True if the ring is a pyran or pyrone ring with a carbonyl group at position 4, False otherwise
    """
    if len(ring_atoms) != 6:
        return False
    
    carbonyl_atom = None
    for atom in ring_atoms:
        if atom.GetAtomicNum() == 8 and len(atom.GetNeighbors()) == 1:
            carbonyl_atom = atom
            break
    
    if not carbonyl_atom:
        return False
    
    neighbor = carbonyl_atom.GetNeighbors()[0]
    if neighbor.GetAtomicNum() != 6 or neighbor.GetIsAromatic():
        return False
    
    ring_atoms.remove(carbonyl_atom)
    for atom in ring_atoms:
        if not (atom.GetAtomicNum() == 6 and atom.GetIsAromatic()):
            return False
    
    return True