"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins – natural pigments containing four pyrrole-like rings
united by four methine bridges to yield an aromatic, conjugated macrocyclic core.
This improved version uses two heuristics:
  1. It searches for a macrocyclic ring between 15 and 18 atoms that contains exactly
     4 nitrogen atoms and in which the vast majority of atoms are aromatic.
  2. It checks if the molecule has a metal center (e.g. Mg, Fe, Zn, Co, Ni, Cu, Pd, Pt)
     that is directly bonded to exactly 4 aromatic nitrogen atoms.
If at least one heuristic is met, the molecule is classified as a porphyrin.
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string using two heuristics:
      Heuristic 1:
        Searches for a macrocyclic ring (15–18 atoms) with exactly 4 nitrogen atoms
        and where at least 80% of the atoms in the ring are aromatic.
      Heuristic 2:
        Checks if the molecule contains a metal atom (e.g. Mg, Fe, Zn, Co, Ni, Cu, Pd, Pt)
        bonded to exactly 4 aromatic nitrogen atoms.
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a porphyrin, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Heuristic 1: Look for a macrocycle (15–18 atoms) with exactly 4 nitrogen atoms.
    # Also require that at least 80% of the atoms in the ring are aromatic.
    for ring in atom_rings:
        ring_size = len(ring)
        if 15 <= ring_size <= 18:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            n_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            aromatic_count = sum(1 for atom in atoms_in_ring if atom.GetIsAromatic())
            # We require exactly 4 nitrogen atoms and that most atoms in the ring are aromatic.
            if n_count == 4 and aromatic_count >= int(0.8 * ring_size):
                reason = (f"Found a macrocyclic ring of {ring_size} atoms with exactly 4 nitrogen atoms "
                          f"and high aromaticity ({aromatic_count}/{ring_size} atoms aromatic), "
                          "consistent with a porphyrin core.")
                return True, reason

    # Heuristic 2: Look for a metal center directly bonded to 4 aromatic nitrogen atoms.
    # List of common metal atomic numbers in porphyrins.
    metal_atomic_nums = {12, 26, 27, 28, 29, 30, 46, 78}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metal_atomic_nums:
            # Count neighboring nitrogen atoms that are aromatic.
            aromatic_n_neighbors = [nbr for nbr in atom.GetNeighbors() 
                                    if nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic()]
            if len(aromatic_n_neighbors) == 4:
                reason = (f"Metal atom {atom.GetSymbol()} (atomic num {atom.GetAtomicNum()}) "
                          "is coordinated to 4 aromatic nitrogen atoms, consistent with a metalloporphyrin core.")
                return True, reason
    
    # If neither heuristic triggers, list both shortcomings.
    reasons = []
    reasons.append("No macrocyclic ring (15–18 atoms with exactly 4 nitrogen atoms and high aromaticity) found.")
    reasons.append("No metal center coordinated to exactly 4 aromatic nitrogen atoms was detected.")
    return False, " ".join(reasons)

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples (e.g. for chlorophyll a')
    test_smiles = ("C1=2N3C(C=C4[N+]5=C(C=C6N7C8=C(C9=[N+](C(=C1)[C@H]([C@@H]9CCC(OC/C=C(/CCC[C@@H](CCC[C@@H](CCCC(C)C)C)C)\C)=O)C)"
                   "[Mg-2]735)[C@@H](C(C8=C6C)=O)C(=O)OC)C(=C4C)C=C)=C(C2C)C=C")
    result, reason = is_porphyrins(test_smiles)
    print("Result:", result)
    print("Reason:", reason)