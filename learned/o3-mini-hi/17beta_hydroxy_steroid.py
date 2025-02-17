"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17β-hydroxy steroid
Definition: A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
The approach is to identify a steroid-like fused tetracyclic system (three six-membered rings and one five-membered ring)
and then look for a hydroxyl group on a saturated (sp3) carbon that is part of one of the 5-membered rings.
Note: Because exact numbering and stereochemistry are difficult to infer automatically, we use the presence of a chiral tag
on the carbon bearing –OH as a proxy for defined beta stereochemistry.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is a 17β-hydroxy steroid based on its SMILES string.
    Our procedure is:
      1. Parse the SMILES.
      2. Obtain ring information and build a connectivity graph of rings (rings are "fused" if they share >= 2 atoms).
      3. Find a connected group of rings (a fused system) that contains at least three 6-membered rings and one 5-membered ring.
      4. Identify a candidate -OH group (on a saturated, non-carbonyl carbon) that is attached on an atom that is part
         of one of the 5-membered rings in this fused system.
      5. Optionally, check that the carbon bearing the –OH has a defined chiral tag.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True and a reason if the molecule appears to be a 17β-hydroxy steroid,
                     False and a reason otherwise.
                     If SMILES parsing fails, returns (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information as a list of tuples, each tuple is atom indices in that ring.
    ring_info = mol.GetRingInfo()
    all_rings = [set(ring) for ring in ring_info.AtomRings()]
    if not all_rings:
        return False, "No rings found in molecule"
    
    # Build a graph of rings: two rings are connected (edge exists) if they share 2 or more atoms.
    n_rings = len(all_rings)
    adjacency = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(all_rings[i].intersection(all_rings[j])) >= 2:
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Find connected components (fused ring systems)
    visited = set()
    components = []
    for i in range(n_rings):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                cur = stack.pop()
                if cur in comp:
                    continue
                comp.add(cur)
                visited.add(cur)
                for neighbor in adjacency[cur]:
                    if neighbor not in comp:
                        stack.append(neighbor)
            components.append(comp)
    
    # Look for a steroid-like fused system:
    steroid_component = None
    for comp in components:
        # Count number of 6-membered and 5-membered rings in the component.
        count6 = 0
        count5 = 0
        for idx in comp:
            ring_size = len(all_rings[idx])
            if ring_size == 6:
                count6 += 1
            elif ring_size == 5:
                count5 += 1
        # For a steroid nucleus we require at least three 6-membered rings and at least one 5-membered ring.
        if count6 >= 3 and count5 >= 1:
            steroid_component = comp
            break
    
    if steroid_component is None:
        return False, "Steroid nucleus (fused ring system with required ring sizes) not found"
    
    # Gather all atom indices in the steroid fused system.
    steroid_atoms = set()
    five_member_rings = []  # also keep a list of five-membered rings in the system
    for idx in steroid_component:
        steroid_atoms.update(all_rings[idx])
        if len(all_rings[idx]) == 5:
            five_member_rings.append(all_rings[idx])
    
    # Define a SMARTS pattern for an aliphatic carbon attached to an -OH (but not part of a carbonyl)
    oh_pattern = Chem.MolFromSmarts("[CX4;!$([C]=O)][OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if not oh_matches:
        return False, "No appropriate hydroxyl (-OH) group on an aliphatic carbon found"
    
    candidate_found = False
    reason = ""
    # Check if any -OH candidate is on an atom that is within our steroid fused system and in a 5-membered ring.
    for match in oh_matches:
        carbon_idx = match[0]  # atom index of the carbon bound to -OH
        if carbon_idx not in steroid_atoms:
            continue  # this OH is not part of the fused steroid nucleus
        # Check if this carbon lies in any of the 5-membered rings of the steroid nucleus
        in_five = any(carbon_idx in ring for ring in five_member_rings)
        if in_five:
            atom = mol.GetAtomWithIdx(carbon_idx)
            if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                candidate_found = True
                reason = ("Steroid nucleus found with a fused system (three 6-membered and one 5-membered rings) "
                          "and a hydroxyl group attached at a chiral carbon in a 5-membered ring (likely C17 with beta-configuration).")
                break
            else:
                candidate_found = True
                reason = ("Steroid nucleus found with a fused system (three 6-membered and one 5-membered rings) "
                          "and a hydroxyl group attached on a carbon in a 5-membered ring (possibly C17), but stereochemistry is not explicitly defined.")
                break

    if not candidate_found:
        return False, "No hydroxyl group found on a carbon in a 5-membered ring of the steroid nucleus (expected for 17β–OH)"
    
    return True, reason