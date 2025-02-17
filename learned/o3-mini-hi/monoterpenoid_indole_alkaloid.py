"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid

A monoterpenoid indole alkaloid is biosynthesised from L-tryptophan and a diisoprenoid
(building block usually derived from secologanin). In our heuristic we require:
  1. The presence of an indole core (free, N-substituted, or via a fused ring system).
  2. The presence of an attached fragment that is not part of the indole core and that
     has a decent number of carbons (roughly 5–15) and is “terpenoid‐like”. A fragment is 
     considered terpenoid‐like if either it contains at least one non‐aromatic double bond or 
     if its carbon atoms are mostly aliphatic (with a high fraction of sp³ hybridization).
  3. Basic overall size limits (total carbon count and molecular weight) are met.
  
Note: This is an approximate heuristic. Some genuine cases might still be missed.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import HybridizationType

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid using a heuristic approach.
    
    Steps:
      1. Attempts to identify an indole core using common SMARTS or by detecting a fused ring system.
      2. For all atoms attached to the indole core, groups them (via a BFS) into connected external fragments.
         Then for each fragment, count the number of carbon atoms, and evaluate its “terpenoid‐likeness”
         by checking if it has a non‐aromatic double bond or a high fraction of sp3 carbons.
      3. Enforce overall size constraints.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple of classification (True/False) and the reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ----- 1. Identify the indole core -----
    # Use SMARTS for free indole and for N-substituted indole.
    indole_free = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c2")
    indole_sub  = Chem.MolFromSmarts("c1ccc2c(c1)[n]c2")
    
    indole_atoms = None
    if mol.HasSubstructMatch(indole_free):
        indole_atoms = set(mol.GetSubstructMatch(indole_free))
    elif mol.HasSubstructMatch(indole_sub):
        indole_atoms = set(mol.GetSubstructMatch(indole_sub))
    else:
        # Try to find a fused ring system that could be an indole
        rings = mol.GetRingInfo().AtomRings()
        for ring1 in rings:
            if len(ring1) == 5:
                # Check if the 5-membered ring has at least one nitrogen
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring1):
                    for ring2 in rings:
                        if len(ring2) == 6:
                            # Check that every atom in the 6-membered ring is aromatic...
                            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring2):
                                # ...and require at least 2 atoms in common for a fused system.
                                if len(set(ring1).intersection(ring2)) >= 2:
                                    indole_atoms = set(ring1).union(ring2)
                                    break
                    if indole_atoms is not None:
                        break
    if indole_atoms is None:
        return False, "No indole core (free, substituted, or fused) found"

    # ----- 2. Look for an attached terpenoid-like fragment -----
    # Find atom indices that are directly attached to the indole core (neighbors in the molecule not in indole)
    external_seeds = set()
    for idx in indole_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in indole_atoms:
                external_seeds.add(nbr_idx)
    
    if not external_seeds:
        return False, "No atoms attached directly to the indole core"

    # Now, get connected components (fragments) of external atoms (using only atoms not in the indole).
    external_fragments = []
    visited = set()
    for seed in external_seeds:
        if seed in visited:
            continue
        # BFS on subgraph of atoms not in indole.
        comp = set()
        queue = [seed]
        while queue:
            current = queue.pop(0)
            if current in comp:
                continue
            comp.add(current)
            atom = mol.GetAtomWithIdx(current)
            for nbr in atom.GetNeighbors():
                nid = nbr.GetIdx()
                if nid not in indole_atoms and nid not in comp:
                    queue.append(nid)
        visited |= comp
        external_fragments.append(comp)
    
    # Evaluate each fragment using two criteria:
    #   - The fragment should contain at least 5 carbon atoms (but not be unreasonably huge, assume <=15 carbons).
    #   - The fragment is assumed "terpenoid-like" if it has at least one non-aromatic double bond
    #     OR if the fraction of sp3 carbons in the fragment is at least 0.4.
    terpene_found = False
    for comp in external_fragments:
        # Count carbon atoms in the fragment.
        comp_carbons = [idx for idx in comp if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        n_carbons = len(comp_carbons)
        if n_carbons < 5 or n_carbons > 15:
            continue

        # Count sp3-hybridized carbon atoms.
        sp3 = 0
        for idx in comp:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and atom.GetHybridization() == HybridizationType.SP3:
                sp3 += 1
        sp3_ratio = sp3 / n_carbons if n_carbons > 0 else 0

        # Now check if any bond (both atoms in the fragment) is a non-aromatic double bond.
        has_nonaromatic_double = False
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            if i in comp and j in comp:
                if bond.GetBondTypeAsDouble() == 2 and (not bond.GetIsAromatic()):
                    has_nonaromatic_double = True
                    break

        # Accept the fragment if either the bond condition holds or if the fragment is mostly sp3.
        if has_nonaromatic_double or (sp3_ratio >= 0.4):
            terpene_found = True
            break

    if not terpene_found:
        return False, "No attached terpene-like fragment (adequate size and/or unsaturation) detected"

    # ----- 3. Overall size and composition checks -----
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, "Not enough carbon atoms to harbor both a tryptophan (indole) and a terpene fragment"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for a typical monoterpenoid indole alkaloid"
    
    return True, "Contains an indole core and an attached terpene-like fragment of adequate size"

# (Optional) Example test calls:
# test_smiles = [
#     "C/C=C\\1/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]",  # Burnamine (should be True)
#     "O=C1N[C@H](C(=O)N[C@H]1CC=2C=3C(NC2C(C=C)(C)C)=CC=C(C3CC=C(C)C)CC=C(C)C)C"               # variecolorin L (False positive previously)
# ]
# for s in test_smiles:
#     print(is_monoterpenoid_indole_alkaloid(s))