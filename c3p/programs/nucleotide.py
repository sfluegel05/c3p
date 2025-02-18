"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is defined as a nucleoside phosphate resulting from the condensation 
of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
This function tests for:
   - Presence of at least one phosphorus atom (phosphate group)
   - A furanose sugar ring (5-membered ring with exactly 1 oxygen and 4 carbons)
   - A nucleobase-like aromatic heterocycle (an aromatic ring with at least 2 nitrogen atoms)
   - The phosphate group attached to the sugar ring via a bridging oxygen atom.
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    A nucleotide should contain a nucleoside (a sugar ring attached to a nucleobase)
    with a phosphate group esterified at the 3' or 5' hydroxy position.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a nucleotide, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of phosphorus atoms (to signal a phosphate group)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found (missing phosphate group)"
    
    # Identify a furanose sugar ring: a 5-membered ring with exactly 1 oxygen and 4 carbons.
    sugar_ring_found = False
    sugar_ring_indices = set()
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:  # looking for 5-membered ring
            num_ox = 0
            num_c = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'O':
                    num_ox += 1
                elif atom.GetSymbol() == 'C':
                    num_c += 1
            if num_ox == 1 and num_c == 4:
                sugar_ring_found = True
                sugar_ring_indices = set(ring)
                break
    
    if not sugar_ring_found:
        return False, "No furanose sugar ring (5-membered ring with 1 oxygen and 4 carbons) found"
    
    # Identify a nucleobase-like aromatic heterocycle:
    # We search for aromatic rings (other than the sugar ring) that have at least 2 nitrogen atoms.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        ring_set = set(ring)
        # Skip the sugar ring (or rings completely included within it)
        if ring_set == sugar_ring_indices or ring_set.issubset(sugar_ring_indices):
            continue
        
        # Accept if the ring is aromatic (all atoms aromatic) and has at least 2 nitrogens
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen >= 2:
                nucleobase_found = True
                break
    
    if not nucleobase_found:
        return False, "No nucleobase-like aromatic heterocycle (ring with ≥2 nitrogens) found"
    
    # Confirm that the phosphate group is attached to the sugar ring.
    # Instead of looking at direct neighbors of phosphorus, we look for sugar ring oxygen atoms
    # that are connected to a phosphorus atom (as in an ester linkage).
    phosphate_attached = False
    for idx in sugar_ring_indices:
        atom = mol.GetAtomWithIdx(idx)
        # Look for oxygen atoms in the sugar ring
        if atom.GetSymbol() == 'O':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 15:  # phosphorus atom found attached to oxygen
                    phosphate_attached = True
                    break
        if phosphate_attached:
            break
    
    if not phosphate_attached:
        return False, "Phosphate group is not attached to the sugar ring via a bridging oxygen atom"
    
    return True, "Molecule contains a nucleoside portion (furanose sugar + nucleobase) with a phosphate group attached"