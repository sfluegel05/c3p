"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogue
A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
This version uses improved heuristics:
  1. The molecule must be a valid RDKit molecule, with molecular weight 80–350 Da.
  2. The molecule must contain at least 8 heavy atoms and ≤35 heavy atoms.
  3. It must have at least 2 nitrogen atoms.
  4. It must be mostly rigid (≤3 rotatable bonds).
  5. The molecule must possess an aromatic heterocyclic core.
     We compute the “aromatic core” as the largest connected set of atoms that are both aromatic and part of a ring.
  6. The aromatic core must be appreciable – it must account for at least 40% of the molecule’s heavy atoms.
  7. We add an extra filter on the aromatic core: its carbon-to-nitrogen ratio should be at least 1.0,
     as typical nucleobases feature a balanced (or carbon‐rich) ring (this helps avoid nitrogen–rich false positives).
     
These rules try to capture both the nearly exclusively aromatic ring system of natural nucleobases and
the presence of substituents seen in nucleobase analogues.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.

    Heuristics used:
      1. Valid RDKit molecule and proper sanitization.
      2. Molecular weight between 80 and 350 Da.
      3. Contains at least 8 and not more than 35 heavy atoms.
      4. Contains at least 2 nitrogen atoms.
      5. Is nearly rigid (≤3 rotatable bonds).
      6. Has an aromatic heterocyclic core: we compute the largest connected set of atoms that are aromatic and in a ring.
      7. The aromatic core must be at least 40% of the heavy atoms.
      8. The aromatic core’s carbon-to-nitrogen ratio must be at least 1.0.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a nucleobase analogue, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"

    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside expected range (80–350 Da)"

    # Count heavy atoms.
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count < 8:
        return False, f"Too few heavy atoms ({heavy_atom_count}); nucleobase analogues are typically larger (≥8 heavy atoms)"
    if heavy_atom_count > 35:
        return False, f"Too many heavy atoms ({heavy_atom_count}); nucleobase analogues are typically small (≤35 heavy atoms)"

    # Count nitrogen atoms (require at least 2).
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogen < 2:
        return False, "Too few nitrogen atoms; a nucleobase analogue should contain at least 2 nitrogen atoms"

    # Count rotatable bonds (molecules should be nearly rigid).
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 3:
        return False, f"Too many rotatable bonds ({n_rotatable}); nucleobase analogues are expected to be rigid (≤3 rotatable bonds)"
    
    # ----- Aromatic core extraction -----
    # We want to find all atoms that are both aromatic and belong to at least one ring.
    ring_info = mol.GetRingInfo()
    aromatic_ring_atoms = set()
    for ring in ring_info.AtomRings():
        # Only consider rings where all atoms are aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_atoms.update(ring)
    if not aromatic_ring_atoms:
        return False, "No aromatic ring atoms detected; nucleobase analogues must be heterocyclic and aromatic"
    
    # Now, among these aromatic ring atoms, we want to group connected atoms into components.
    # We build a graph between aromatic atoms based on bonds between them.
    aromatic_atoms = list(aromatic_ring_atoms)
    # Create a mapping for connectivity:
    neighbors = {idx: set() for idx in aromatic_atoms}
    for idx in aromatic_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in aromatic_ring_atoms:
                neighbors[idx].add(nbr_idx)
    
    # Perform DFS to get connected components.
    def dfs(start, seen):
        stack = [start]
        comp = set()
        while stack:
            node = stack.pop()
            if node in comp:
                continue
            comp.add(node)
            for nbr in neighbors[node]:
                if nbr not in comp:
                    stack.append(nbr)
        return comp

    components = []
    visited = set()
    for idx in aromatic_atoms:
        if idx not in visited:
            comp = dfs(idx, visited)
            visited.update(comp)
            components.append(comp)
    
    # Identify the largest aromatic component.
    largest_core = max(components, key=len)
    core_size = len(largest_core)
    aromatic_core_fraction = core_size / heavy_atom_count
    if aromatic_core_fraction < 0.4:
        return False, ("The largest aromatic core contains only {} out of {} heavy atoms (fraction {:.2f}); "
                       "nucleobase analogues are expected to be predominantly aromatic".format(core_size, heavy_atom_count, aromatic_core_fraction))
    
    # Evaluate the carbon-to-nitrogen ratio in the aromatic core.
    core_n_C = 0
    core_n_N = 0
    for idx in largest_core:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            core_n_C += 1
        elif atom.GetAtomicNum() == 7:
            core_n_N += 1
    # To avoid division by zero, if there are no nitrogens in the core, fail.
    if core_n_N == 0:
        return False, "Aromatic core contains no nitrogen atoms; nucleobase analogues must be heterocyclic"
    core_c_to_n = core_n_C / core_n_N
    if core_c_to_n < 1.0:
        return False, (f"Aromatic core C/N ratio of {core_c_to_n:.2f} is too low; nucleobase analogues typically are not excessively nitrogen‐rich")
    
    return True, ("Molecule exhibits features of a nucleobase analogue: appropriate molecular weight, nitrogen content, "
                  "small size, rigidity, and a predominantly aromatic heterocyclic core with balanced C/N ratio.")

# Example testing section (uncomment to test several cases)
# if __name__ == "__main__":
#     test_cases = [
#         # Expected True examples (nucleobase analogues)
#         ("C=12C(=NC=NC1NCC(=O)N)NC=N2", "N(6)-carbamoylmethyladenine"),
#         ("Nc1nc(=O)[nH]cc1CO", "5-(hydroxymethyl)cytosine"),
#         ("Nc1ncc2nc[nH]c2n1", "2-aminopurine"),
#         ("O=c1[nH]cnc2nc[nH]c12", "hypoxanthine"),
#         ("C=12C(=NC=NC1N)NC(N2)=O", "8-oxoadenine"),
#         # Expected False examples (non-nucleobase analogues)
#         ("c1ncnnn1", "1,2,3,5-tetrazine"),
#         ("C1=CN=C(C=N1)Cl", "2-chloropyrazine"),
#         ("O(C1=NC=CN=C1)C", "Methoxypyrazine"),
#     ]
#     for smi, name in test_cases:
#         res, reason = is_nucleobase_analogue(smi)
#         print(f"SMILES: {smi}\nName: {name}\nResult: {res}\nReason: {reason}\n")