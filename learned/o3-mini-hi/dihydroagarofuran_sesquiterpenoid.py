"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: Dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.
The dihydroagarofuran skeleton is expected to be a fully saturated (non‐aromatic)
tricyclic core built from a fused 5–6–5 ring system that, after removal of decorating substituents
(using the Murcko scaffold), contains roughly 15 carbon atoms (allowing 14–16, or one less if one oxygen is present)
and at most one oxygen.
"""
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors
from itertools import combinations

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    Approach:
      1. Parse the SMILES string.
      2. Compute the Murcko scaffold to remove peripheral substituents.
      3. Use GetSymmSSSR to obtain all ring systems in the scaffold.
      4. Filter rings: only those of size 5 or 6 and entirely non‐aromatic.
      5. Consider all combinations of three candidate rings.
         For each triplet:
           - Check whether at least two pairs share at least 2 atoms (fused by two-atom bonds).
           - Check that the sorted ring sizes equal [5, 5, 6] regardless of order.
           - Compute the union of atoms in the triplet; from this union, count only carbons and oxygens.
           - For a valid dihydroagarofuran, if there is no oxygen then the number of carbons should be between 14 and 16,
             and if there is one oxygen then between 13 and 15.
    
    Returns:
       bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise.
       str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Compute the Murcko scaffold.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not compute Murcko scaffold."
    
    # Compute ring systems using GetSymmSSSR (more robust than GetRingInfo at times).
    ssr = Chem.GetSymmSSSR(scaffold)
    if not ssr:
        return False, "No rings detected in the Murcko scaffold."
    
    candidate_rings = []
    # Loop over each ring (ssr returns a tuple of atom index tuples)
    for ring in ssr:
        ring_atoms = list(ring)
        # Only consider rings of size 5 or 6
        if len(ring_atoms) not in (5, 6):
            continue
        # Check that every atom in the ring is non-aromatic.
        if any(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
            continue
        candidate_rings.append(ring_atoms)
    
    if len(candidate_rings) < 3:
        return False, (f"Only {len(candidate_rings)} candidate ring(s) (of size 5 or 6 and non‐aromatic) were detected in the scaffold; "
                       "at least 3 are needed for a fused 5–6–5 core.")
    
    # Consider each combination of 3 candidate rings.
    for triplet in combinations(candidate_rings, 3):
        fused_edges = 0
        # Count how many ring pairs share at least 2 atoms.
        for ring_a, ring_b in combinations(triplet, 2):
            if len(set(ring_a).intersection(ring_b)) >= 2:
                fused_edges += 1
        # For a proper fused system, we require at least 2 fused edges.
        if fused_edges < 2:
            continue
        
        # Check that the ring sizes (order independent) are exactly two five-membered and one six-membered ring.
        ring_sizes = sorted([len(r) for r in triplet])
        if ring_sizes != [5, 5, 6]:
            continue
        
        # Get the union of atom indices in these three rings.
        union_atoms = set()
        for ring in triplet:
            union_atoms.update(ring)
        
        # Restrict attention to atoms that are carbon (atomic number 6) or oxygen (atomic number 8) 
        # (ignoring any heteroatoms that might have survived).
        core_atoms = []
        valid_core = True
        for idx in union_atoms:
            atom = scaffold.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() in (6, 8):
                core_atoms.append(atom)
            else:
                valid_core = False
                break
        if not valid_core:
            continue
        
        # Count carbons and oxygens.
        c_count = sum(1 for atom in core_atoms if atom.GetAtomicNum() == 6)
        o_count = sum(1 for atom in core_atoms if atom.GetAtomicNum() == 8)
        
        if o_count > 1:
            continue
        
        # Check carbon count thresholds.
        if o_count == 0:
            if not (14 <= c_count <= 16):
                continue
        else:  # o_count == 1
            if not (13 <= c_count <= 15):
                continue
        
        # Found a matching fused 5-6-5 core.
        return True, ("Molecule's scaffold contains a fused 5–6–5 tricyclic non‐aromatic core with "
                      f"{c_count} carbon(s) and {o_count} oxygen(s), consistent with a dihydroagarofuran skeleton.")
    
    return False, "No fused 5–6–5 tricyclic core (non‐aromatic, with ~15 carbons and at most 1 oxygen) was detected in the molecule's scaffold."

# For testing purposes, you might call the function as follows:
if __name__ == "__main__":
    # Example: wilfordinine B SMILES
    smiles_example = "[C@@]1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(C)=O)[C@@]([H])([C@H]([C@H]([C@@]3([C@@H](OC(=O)C)[C@H]2OC(C)=O)COC(=O)C)OC(=O)C)OC(C)=O)[C@@](COC(C5=C([C@H]1C)C=CN=C5)=O)(C)O4)C)[H])=O)(C)O"
    result, reason = is_dihydroagarofuran_sesquiterpenoid(smiles_example)
    print("Result:", result)
    print("Reason:", reason)