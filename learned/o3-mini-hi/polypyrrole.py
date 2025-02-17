"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
Definition: A compound composed of two or more pyrrole units.
A pyrrole unit is defined as a five-membered aromatic ring containing exactly one nitrogen.
"""

from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a molecule that contains two or more pyrrole units.
    
    A pyrrole unit is defined as a five-membered aromatic ring containing exactly one nitrogen atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a polypyrrole, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Counter for pyrrole rings detected.
    pyrrole_count = 0
    
    # Iterate over each ring in the molecule.
    for ring in atom_rings:
        # Consider only rings of exactly 5 atoms.
        if len(ring) != 5:
            continue
        
        # Get atoms in the ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # All atoms in the ring must be aromatic.
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        
        # Count nitrogen atoms in the ring.
        n_nitrogen = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
        
        # For a pyrrole unit, exactly one nitrogen is expected in a 5-membered aromatic ring.
        if n_nitrogen == 1:
            pyrrole_count += 1
    
    if pyrrole_count >= 2:
        return True, f"Contains {pyrrole_count} pyrrole units (5-membered aromatic rings with one nitrogen)"
    else:
        return False, f"Contains {pyrrole_count} pyrrole unit(s); need at least two for a polypyrrole"
        
# Example test cases:
if __name__ == "__main__":
    # These are example SMILES strings; one could test with known polypyrroles here.
    test_smiles = [
        "C1=CC2=C(C=C1)C=CC=C2",  # Not polypyrrole (benzene rings)
        "c1cc[nH]c1",            # Single pyrrole unit
        "c1cc[nH]c1-c2cc[nH]c2"   # Two pyrrole units connected
    ]
    
    for s in test_smiles:
        result, reason = is_polypyrrole(s)
        print(f"SMILES: {s}\nResult: {result} \nReason: {reason}\n")