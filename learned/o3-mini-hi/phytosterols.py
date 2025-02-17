"""
Classifies: CHEBI:26125 phytosterols
"""
#!/usr/bin/env python
"""
Classifies: Phytosterols – Sterols similar to cholesterol that occur in plants.
Phytosterols are characterized by a steroid nucleus (a fused system of three six–membered rings 
and one five–membered ring) and typically a hydroxyl group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import itertools

def is_phytosterols(smiles: str):
    """
    Determines if a molecule belongs to the phytosterols class based on its SMILES string.
    
    Phytosterols are plant sterols defined as sterols having a fused four–ring (steroid)
    system (three six–membered rings and one five–membered ring) and usually a hydroxyl group,
    typically at the 3β position. In this implementation we:
      1. Parse the SMILES string.
      2. Retrieve ring information and consider rings that are either 5 or 6 members.
      3. Look for a combination of 4 rings that together yield a fused network with about 16–18 unique carbon atoms,
         with exactly one five–membered ring.
      4. Check for the presence of at least one hydroxyl group.
      5. Verify the molecular weight is in the expected range for sterols (~300–600 Da).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phytosterol, False otherwise.
        str: Explanation (reason) for the classification decision.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information without checking for IsInitialized (since it is always available)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is represented as a tuple of atom indices

    # Filter rings that are 5 or 6 members (typical sizes for the steroid nucleus)
    candidate_rings = [r for r in rings if len(r) in (5, 6)]
    if len(candidate_rings) < 4:
        return False, "Not enough candidate rings (size 5 or 6) to form a steroid nucleus"

    # Look for a combination of 4 rings that may form the steroid nucleus.
    # We require:
    #  - The united rings contain about 16-18 unique carbon atoms.
    #  - Exactly one of the selected rings is a five-membered ring.
    #  - The rings are fused (for at least one pair, the intersection has at least 2 atoms).
    nucleus_found = False
    for combo in itertools.combinations(candidate_rings, 4):
        unique_atoms = set().union(*combo)
        # Count only carbon atoms in the union of these rings
        carbon_count = sum(1 for idx in unique_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if not (16 <= carbon_count <= 18):
            continue
        five_membered = sum(1 for ring in combo if len(ring) == 5)
        if five_membered != 1:
            continue
        # Build a simple connectivity graph among the four rings.
        # Two rings are considered fused if they share at least 2 atoms.
        graph = {i: set() for i in range(4)}
        for i, j in itertools.combinations(range(4), 2):
            if len(set(combo[i]).intersection(combo[j])) >= 2:
                graph[i].add(j)
                graph[j].add(i)
        # Check if the graph is connected (all rings belong to one fused network)
        visited = set()
        nodes = [0]
        while nodes:
            current = nodes.pop()
            if current in visited:
                continue
            visited.add(current)
            nodes.extend(graph[current] - visited)
        if len(visited) == 4:
            nucleus_found = True
            break

    if not nucleus_found:
        return False, "Molecule does not contain a characteristic fused steroid nucleus (4 rings with 3 6-membered and 1 5-membered fused rings with ~17 carbons)"

    # Check for presence of at least one hydroxyl group (–OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxyl_pattern is None:
        return False, "Error creating hydroxyl SMARTS pattern"
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Molecule lacks a hydroxyl group (–OH), a hallmark of sterols"

    # Verify the molecular weight is within a typical range for sterols (approximately 300-600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the typical range for phytosterols"
    
    return True, "Molecule contains a fused steroid nucleus with ~17 carbon atoms and a hydroxyl group typical of phytosterols"


# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with one of the provided examples: ergosta-5,7-dien-3beta-ol
    test_smiles = "[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@H](C)C(C)C"
    result, reason = is_phytosterols(test_smiles)
    print(result, reason)