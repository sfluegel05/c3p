"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid

A monoterpenoid indole alkaloid is biosynthesised from L-tryptophan and a diisoprenoid
(building block usually derived from secologanin). Our improved heuristic requires:
  1. Presence of an indole core (free, N-substituted, or fused).
  2. An attached external fragment (via the indole core) that has a reasonable number 
     of carbons (roughly 5–20) and is terpenoid-like: it should have at least one non‐aromatic 
     double bond or be largely aliphatic (with a high fraction of sp3 carbons), a high carbon
     fraction among heavy atoms (>0.7) and low aromatic carbon fraction (<0.3).
  3. Overall size limits (sufficient total carbon count and molecular weight).
  
Note: This is an approximate heuristic and may miss or misclassify borderline cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import HybridizationType

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid using an improved heuristic.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): Tuple classification and reason.
    """
    # ----- 0. Parse SMILES and basic error check -----
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ----- 1. Identify the indole core -----
    # Use two SMARTS: free indole and N-substituted indole.
    indole_free = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c2")
    indole_sub  = Chem.MolFromSmarts("c1ccc2c(c1)[n]c2")
    
    indole_atoms = None
    if mol.HasSubstructMatch(indole_free):
        indole_atoms = set(mol.GetSubstructMatch(indole_free))
    elif mol.HasSubstructMatch(indole_sub):
        indole_atoms = set(mol.GetSubstructMatch(indole_sub))
    else:
        # As backup: look for fused ring systems that could form an indole
        rings = mol.GetRingInfo().AtomRings()
        for ring1 in rings:
            if len(ring1) == 5 and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring1):
                for ring2 in rings:
                    if len(ring2) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring2):
                        if len(set(ring1).intersection(ring2)) >= 2:
                            indole_atoms = set(ring1).union(ring2)
                            break
                if indole_atoms is not None:
                    break
    if indole_atoms is None:
        return False, "No indole core (free, substituted, or fused) found"

    # ----- 2. Identify external attached fragments that could be the terpene-like part -----
    # First: determine atom indices directly attached to the indole core
    external_seeds = set()
    for idx in indole_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in indole_atoms:
                external_seeds.add(nbr_idx)
    
    if not external_seeds:
        return False, "No atoms directly attached to the indole core"

    # Now group external atoms into connected fragments (only atoms NOT in the indole core)
    external_fragments = []
    visited = set()
    for seed in external_seeds:
        if seed in visited:
            continue
        comp = set()
        queue = [seed]
        while queue:
            cur = queue.pop(0)
            if cur in comp:
                continue
            comp.add(cur)
            atom = mol.GetAtomWithIdx(cur)
            for nbr in atom.GetNeighbors():
                nid = nbr.GetIdx()
                if nid not in indole_atoms and nid not in comp:
                    queue.append(nid)
        visited |= comp
        external_fragments.append(comp)
    
    # Now evaluate each external fragment for terpenoid criteria
    terpene_found = False
    for frag in external_fragments:
        # Count carbon atoms in the fragment
        frag_carbons = [idx for idx in frag if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        n_carbons = len(frag_carbons)
        # Accept fragments between 5 and 20 carbons (wider range to catch borderline cases)
        if n_carbons < 5 or n_carbons > 20:
            continue
        
        # Calculate fraction of sp3-hybridized carbons
        sp3 = 0
        aromatic_carbons = 0
        total_heavy = 0  # count non-H atoms in fragment
        for idx in frag:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 1:  # ignore hydrogens
                total_heavy += 1
            if atom.GetAtomicNum() == 6:
                if atom.GetHybridization() == HybridizationType.SP3:
                    sp3 += 1
                if atom.GetIsAromatic():
                    aromatic_carbons += 1
        sp3_ratio = sp3 / n_carbons if n_carbons > 0 else 0
        aromatic_ratio = aromatic_carbons / n_carbons if n_carbons > 0 else 0

        # Check heavy atom composition: fraction of carbons among non-hydrogens
        frag_carbon_count = sum(1 for idx in frag if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        frag_heavy_count = sum(1 for idx in frag if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
        carbon_fraction = frag_carbon_count / frag_heavy_count if frag_heavy_count > 0 else 0

        # Look for at least one non‐aromatic double bond within the fragment.
        has_nonaromatic_double = False
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            if i in frag and j in frag:
                # Check that the bond is double and not aromatic.
                if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                    has_nonaromatic_double = True
                    break

        # We expect the external terpene fragment to be largely aliphatic.
        # Thus we require a decent sp3 or the presence of a non‐aromatic double bond,
        # a high fraction of carbons among heavy atoms (>0.7) and low aromaticity (<0.3).
        if (has_nonaromatic_double or sp3_ratio >= 0.4) and carbon_fraction >= 0.7 and aromatic_ratio < 0.3:
            terpene_found = True
            break

    if not terpene_found:
        return False, "No attached terpene-like fragment (adequate size, unsaturation, and aliphatic character) detected"

    # ----- 3. Overall size and composition checks -----
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, "Not enough carbon atoms to support both an indole core and a terpene-like fragment"

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for a typical monoterpenoid indole alkaloid"
    
    return True, "Contains an indole core and an attached terpene-like fragment of adequate size and aliphatic character"

# (Optional) Test example calls:
# test_smiles = [
#     "COC([C@@H]1C[C@@]23CCC[N+]4(CCC5(C6=CC=CC=C6N(C(=O)OC)[C@]15CC2)[C@]34[H])[O-])=O",  # Pleiocarpoline, expected True
#     "C/C=C\\1/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]"   # Burnamine, expected True
# ]
# for s in test_smiles:
#     print(is_monoterpenoid_indole_alkaloid(s))