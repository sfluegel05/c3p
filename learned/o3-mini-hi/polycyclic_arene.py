"""
Classifies: CHEBI:33848 polycyclic arene
"""
#!/usr/bin/env python
"""
Classifies: Polycyclic Arene, defined as a polycyclic aromatic hydrocarbon.
A polycyclic arene must contain at least two fused candidate rings.
Candidate rings are defined as rings in which every atom is aromatic and is a carbon.
In addition, since polycyclic aromatic hydrocarbons (PAHs) are almost entirely made of carbons and hydrogens,
we require that the overall heavy-atom composition (ignoring hydrogens) is strongly dominated by carbon.
This extra criterion helps filter false positives arising from molecules that have only a small fused aromatic
substructure amidst many heteroatoms.
"""

from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (a PAH)
    based on its SMILES string.
    
    The algorithm works by:
      1. Converting the SMILES string to an RDKit molecule.
      2. Extracting ring information and flagging rings as candidate rings:
         a candidate ring is one where every atom is aromatic and is carbon.
      3. Building a connectivity graph among rings (two rings are fused if they share at least 2 atoms).
      4. Finding connected components among rings and flagging any component that contains
         at least two candidate rings.
      5. Finally, we check that the molecule’s heavy atoms (non‐hydrogen atoms) are dominated by carbon.
         A ratio of carbon atoms to heavy atoms of at least 0.85 is required.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as a
                      polycyclic aromatic hydrocarbon. The second element is a string giving
                      the reason.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring is represented as a tuple of atom indices)
    ring_info = mol.GetRingInfo()
    atom_rings = [set(ring) for ring in ring_info.AtomRings()]
    if not atom_rings:
        return False, "No rings detected in the molecule"
    
    # Define helper to decide if a ring is an aromatic carbon ring.
    def is_aromatic_carbon_ring(ring):
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check that the atom is aromatic and is carbon (atomic number 6)
            if not atom.GetIsAromatic() or atom.GetAtomicNum() != 6:
                return False
        return True
    
    # Flag each ring if it is a candidate (fully aromatic and all carbon)
    candidate_flags = [is_aromatic_carbon_ring(ring) for ring in atom_rings]
    
    # Build a ring connectivity graph.
    # Two rings are fused if they share at least 2 atoms.
    n = len(atom_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if len(atom_rings[i].intersection(atom_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components in the ring graph.
    visited = set()
    for i in range(n):
        if i in visited:
            continue
        # Depth-first search to get the full connected component.
        stack = [i]
        component = set()
        while stack:
            node = stack.pop()
            if node not in component:
                component.add(node)
                stack.extend(ring_graph[node] - component)
        visited |= component
        
        # Count candidate rings in this component.
        candidate_count = sum(1 for idx in component if candidate_flags[idx])
        if candidate_count >= 2:
            # We now perform an additional check on the overall composition.
            # PAHs are mainly hydrocarbons.
            heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]  # exclude H's
            if not heavy_atoms:
                return False, "No heavy atoms found"
            total_heavy = len(heavy_atoms)
            carbon_count = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
            carbon_ratio = carbon_count / total_heavy
            
            if carbon_ratio < 0.85:
                return False, ("Found a fused aromatic candidate system but the molecule's heavy-atom composition " +
                               f"is only {carbon_ratio:.2f} carbon, suggesting significant heteroatom content")
            # Optionally, one might also require that the fused candidate system covers a large fraction of the molecule.
            # For example, one may compute the union of candidate ring atoms and compare to the total heavy atoms.
            component_candidate_atoms = set()
            for idx in component:
                if candidate_flags[idx]:
                    component_candidate_atoms |= atom_rings[idx]
            if len(component_candidate_atoms) < 10:
                # Many simple substituent aromatic rings can have only 6-8 atoms total; we choose 10 as a heuristic.
                return False, ("Fused candidate system found but overall number of ring atoms is low " +
                               f"({len(component_candidate_atoms)}), so the molecule likely is not a large PAH")
            
            return True, f"Found a fused aromatic hydrocarbon system with {candidate_count} candidate rings"
    
    return False, "No fused aromatic hydrocarbon system detected"


# Example testing code (you may remove this block when integrating into a larger system)
if __name__ == '__main__':
    test_smiles = [
        "c1ccc2c(c1)ccc1cc3c(ccc4ccccc34)cc21",  # dibenz[a,h]anthracene (should be True)
        "c1ccc2cc3ccccc3cc2c1",                  # anthracene (should be True)
        "Oc1cc2ccc3cc4ccccc4c4ccc(c1)c2c34",      # 2-Hydroxybenzo[a]pyrene (should be True)
        "O1C(OC=2C3=C(OC)C=C(OC)C=C3C=C4C2CO[C@H](C4)C)C(O)C(O)C(C1CO)OC",  # An example false positive
    ]
    for smi in test_smiles:
        result, reason = is_polycyclic_arene(smi)
        print(smi, "=>", result, ":", reason)