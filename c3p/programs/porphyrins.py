"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins – natural pigments that contain a macrocyclic core
built from four pyrrole-like rings united by four methine bridges.
This program uses two heuristics:
  1. It checks for a large ring (15–18 atoms) that has exactly 4 nitrogen atoms.
  2. It checks if the molecule contains a metal atom (e.g. Mg, Fe, Zn, Co, Pd, Pt)
     that is coordinated to exactly 4 nitrogen atoms.
If either criterion is met, the molecule is classified as a porphyrin.
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Heuristic 1:
      Searches the ring information for a large (approx. 15–18 atoms) ring that has exactly 4 nitrogen atoms.
    Heuristic 2:
      Checks if the molecule contains a metal that is bonded to exactly 4 nitrogen atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a porphyrin, False otherwise.
        str: Reason explaining the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Heuristic 1: Look for a macrocycle (15–18 atoms) with exactly 4 nitrogens.
    macrocycle_found = False
    macrocycle_reason = ""
    for ring in atom_rings:
        ring_size = len(ring)
        # Look for a ring in the range 15 to 18 atoms (porphyrin macrocycle is typically 16 atoms)
        if 15 <= ring_size <= 18:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count == 4:
                macrocycle_found = True
                macrocycle_reason = (f"Found a macrocyclic ring of {ring_size} atoms "
                                     f"with exactly {n_count} nitrogen atoms, consistent with a porphyrin core.")
                break
    if macrocycle_found:
        return True, macrocycle_reason

    # Heuristic 2: Look for a metal center coordinated to 4 nitrogen atoms.
    # List of common metal atomic numbers found in porphyrins:
    # Magnesium (12), Iron (26), Cobalt (27), Nickel (28), Copper (29), Zinc (30), Palladium (46), Platinum (78)
    metal_atomic_nums = {12, 26, 27, 28, 29, 30, 46, 78}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metal_atomic_nums:
            # Count neighboring nitrogen atoms (atom must be directly bonded).
            n_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 7]
            if len(n_neighbors) == 4:
                return True, (f"Metal atom {atom.GetSymbol()} (atomic num {atom.GetAtomicNum()}) "
                              "is coordinated to 4 nitrogen atoms, consistent with a metalloporphyrin core.")
    
    # If neither heuristic triggers, we do not classify as a porphyrin.
    reasons = []
    if not macrocycle_found:
        reasons.append("No macrocyclic ring (15–18 atoms with exactly 4 nitrogen atoms) found.")
    reasons.append("No metal center coordinated to exactly 4 nitrogen atoms was detected.")
    return False, " ".join(reasons)

# Example usage:
if __name__ == "__main__":
    # Example: using one of the provided SMILES, e.g. for chlorophyll a'
    test_smiles = ("C1=2N3C(C=C4[N+]5=C(C=C6N7C8=C(C9=[N+](C(=C1)[C@H]([C@@H]9CCC(OC/C=C(/CCC[C@@H](CCC[C@@H](CCCC(C)C)C)C)\C)=O)C)"
                   "[Mg-2]735)[C@@H](C(C8=C6C)=O)C(=O)OC)C(=C4C)C=C)=C(C2C)C=C")
    result, reason = is_porphyrins(test_smiles)
    print("Result:", result)
    print("Reason:", reason)