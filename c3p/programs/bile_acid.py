"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Hydroxy-5β-cholanic acids (bile acids)
Definition: Bile acids are hydroxy-5β-cholanic acids (or their glycine/taurine conjugates). They typically
display a fused steroid nucleus (at least one five-membered ring plus three six-membered rings with significant overlap),
a carboxyl (or conjugated amide) group on an aliphatic side-chain – that is, the acid carbon must be located away 
(from the steroid core by at least 3 bonds) – and extra hydroxyl substituents.
The algorithm below applies SMARTS patterns for a free acid, glycine conjugate, or taurine conjugate, counts hydroxyls,
checks approximate molecular weight (350–700 Da), and then determines if a fused steroid core is present.
Finally, the distance (shortest path) from the acid carbon to any core atom must be at least 3 bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Classify whether a molecule is a bile acid (hydroxy-5β-cholanic acid or its glycine/taurine conjugate).
    
    The molecule must pass these tests:
      1. It is parsed correctly.
      2. It must contain either a free acid group OR an amide bond found in glycine or taurine conjugates.
         Three SMARTS patterns are used.
      3. It must show at least three hydroxyl groups ([OX2H]).
      4. Its molecular weight should be roughly between 350 and 700 Da.
      5. It must have a fused steroid nucleus as defined by at least one 5-membered ring and three 6-membered rings.
         Moreover, when collecting all atoms in any 5- or 6-membered ring, the fused core contains at least 15 atoms.
      6. The acid (or conjugate acid) group must be on a side-chain: the acid carbon (the carbonyl carbon)
         must not be part of the fused steroid core and the shortest path from the acid carbon to the steroid core
         must be at least 3 bonds.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if classified as a bile acid; False otherwise.
        str: Reason for the decision.
    """
    
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the acid/conjugation group.
    # Free carboxylic acid (C(=O)[O-] or neutral –OH acid) not in a ring.
    acid_pattern = Chem.MolFromSmarts("[$([CX3](=O)[O;H,-]);!R]")
    # Glycine conjugate: an amide where the acyl carbon is attached to –NCC(=O)[O;H,-]
    glycine_pattern = Chem.MolFromSmarts("[$([CX3](=O)NCC(=O)[O;H,-]);!R]")
    # Taurine conjugate: an amide where the acyl carbon is attached to –NCCS(=O)(=O)[O;H,-]
    taurine_pattern = Chem.MolFromSmarts("[$([CX3](=O)NCCS(=O)(=O)[O;H,-]);!R]")
    
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    taurine_matches = mol.GetSubstructMatches(taurine_pattern)
    
    if not (acid_matches or glycine_matches or taurine_matches):
        return False, "Missing free acid or recognized conjugated acid group essential for cholanic acids"
    
    # Check for hydroxyl groups ([OX2H]). Require at least 3 hydroxyls.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Too few hydroxyl substituents (found {len(hydroxyl_matches)}; expect at least 3)"
    
    # Check molecular weight in the rough range of bile acids (350-700 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 700):
        return False, f"Molecular weight {mol_wt:.1f} out of expected bile acid range (350-700 Da)"
    
    # Analyze ring systems to determine if a suitable fused steroid nucleus is present.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found; expected a steroid nucleus"
    
    # Count the rings: bile acids generally have at least one 5-membered ring and three 6-membered rings.
    count_5 = sum(1 for ring in rings if len(ring) == 5)
    count_6 = sum(1 for ring in rings if len(ring) == 6)
    if count_5 < 1 or count_6 < 3:
        return False, f"Ring system does not meet steroid criteria (found {count_5} 5-membered and {count_6} 6-membered rings)"
    
    # Assemble all atoms that are part of a 5- or 6-membered ring.
    steroid_core_atoms = set()
    for ring in rings:
        if len(ring) in (5, 6):
            steroid_core_atoms.update(ring)
    if len(steroid_core_atoms) < 15:
        return False, "Fused steroid core not clearly discernible (not enough atoms in 5-/6-membered rings)"
    
    # Check that the acid (or conjugate) group is on a side-chain away from the steroid core.
    # Our test: for any acid match, the acid carbon (first atom in the SMARTS match) must NOT be in the steroid core,
    # and its shortest path distance to any atom in the steroid core must be at least 3 bonds.
    
    # Combine all acid-like matches
    all_acid_matches = list(acid_matches) + list(glycine_matches) + list(taurine_matches)
    
    sidechain_ok = False
    for match in all_acid_matches:
        acid_carbon = match[0]  # first atom in the match; expected to be the acyl (C=O) carbon.
        if acid_carbon in steroid_core_atoms:
            continue  # Acid carbon is part of the core; not acceptable.
        # Compute the minimum shortest-path distance from the acid carbon to any atom in the steroid core.
        min_distance = float("inf")
        for core_atom in steroid_core_atoms:
            sp = Chem.GetShortestPath(mol, acid_carbon, core_atom)
            if sp and (len(sp) - 1) < min_distance:
                min_distance = len(sp) - 1
        # Adjusted threshold: require at least 3 bonds separation.
        if min_distance >= 3:
            sidechain_ok = True
            break

    if not sidechain_ok:
        return False, "Acid/conjugate group appears too close to the steroid core; expected a side chain of at least 3 bonds"
    
    return True, "Molecule shows a fused steroid nucleus with a separated acid (or conjugated amide) side chain and extra hydroxyls typical for bile acids"

# For testing (this block is optional):
if __name__ == "__main__":
    test_smiles = "O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(C[C@H](O)CC4)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CCC(O)=O)C)[H])[H])C"
    result, reason = is_bile_acid(test_smiles)
    print(result, reason)