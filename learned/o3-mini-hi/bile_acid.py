"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Hydroxy-5β-cholanic acids (bile acids)
Definition: Any member of a group of hydroxy-5β-cholanic acids occurring in bile, where they are present 
as the sodium salts of their amides with glycine or taurine. In mammals, bile acids almost invariably have 
the 5β-configuration. Because the 5β stereochemistry is hard to verify by a generic SMARTS search, here we 
approximate a bile acid by requiring: 
  • a steroid nucleus (approximated as a fused ring system with one 5-membered and at least three 6-membered rings),
  • an acid function (either a free carboxylic acid or an amide expected in glycine/taurine conjugates),
  • an “acid side chain” (the acid or amide group must be separated from the fused poly‐ring system by at least a few bonds),
  • and at least one additional hydroxyl group.
Note: This function uses several SMARTS patterns and distance‐based checks to improve specificity.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid (hydroxy-5β-cholanic acid or its glycine/taurine conjugate)
    based on its SMILES string.
    
    To be classified as a bile acid the molecule must pass several tests:
      1. It must be parsed correctly.
      2. It must contain either a free carboxylic acid group OR show an amide bond arising from conjugation
         with glycine or taurine. (We check for three patterns below.)
      3. It must have at least one extra hydroxyl substituent (apart from the acid/conjugate OH).
      4. Its molecular weight should lie roughly in the 350–700 Da range.
      5. It must show a fused steroid nucleus. Here we require that at least one five‐membered ring and three
         six‐membered rings (from the set of all rings) are present.
      6. Finally, we “localize” the acid (or its conjugated amide) to a side‐chain: the attachment point (the carbon 
         of the acid group that binds the chain) should be non‐ring (or at least not directly fused into the steroid core).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a bile acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS patterns:
    # Free carboxylic acid: matches –C(=O)[O;H,-]
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    # Glycine conjugate: the acid is amidated with a glycine fragment: –C(=O)NCC(=O)[O;H,-]
    glycine_pattern = Chem.MolFromSmarts("C(=O)NCC(=O)[O;H,-]")
    # Taurine conjugate: the acid is amidated with a taurine fragment: –C(=O)NCCS(=O)(=O)[O;H,-]
    taurine_pattern = Chem.MolFromSmarts("C(=O)NCCS(=O)(=O)[O;H,-]")

    # Must have at least one of the acid (free or conjugate) patterns.
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    taurine_matches = mol.GetSubstructMatches(taurine_pattern)
    if not (acid_matches or glycine_matches or taurine_matches):
        return False, "Missing carboxylic acid or conjugated acid group essential for cholanic acids"
    
    # Check for hydroxyl groups ([OX2H]). We require at least 2 – one may belong to the acid motif.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Too few hydroxyl substituents (expect at least one extra hydroxyl besides acid/conjugate OH)"

    # Check molecular weight (rough guideline for bile acids is 350-700 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 700):
        return False, f"Molecular weight {mol_wt:.1f} out of expected bile acid range"
    
    # Check for a fused steroid nucleus:
    # We require that the molecule has rings and that among its rings there is (approximately) one 5-membered and three 6-membered rings.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found, expected a steroid nucleus"
    
    count_5 = 0
    count_6 = 0
    for ring in rings:
        if len(ring) == 5:
            count_5 += 1
        elif len(ring) == 6:
            count_6 += 1
    if count_5 < 1 or count_6 < 3:
        return False, f"Ring system does not match steroid nucleus criteria (found {count_5} 5-membered and {count_6} 6-membered rings)"
    
    # For further specificity, collect atoms that are part of any 5- or 6-membered ring; assume these form the steroid core.
    steroid_core_atoms = set()
    for ring in rings:
        if len(ring) in (5, 6):
            for atom_idx in ring:
                steroid_core_atoms.add(atom_idx)
    if len(steroid_core_atoms) < 15:
        return False, "Fused steroid core not clearly discernible"
    
    # Now check that the acid (or conjugate) group is on a side-chain beyond the steroid nucleus.
    # We combine all matches from acid, glycine and taurine patterns.
    combined_matches = list(acid_matches) + list(glycine_matches) + list(taurine_matches)
    sidechain_ok = False
    for match in combined_matches:
        # In our SMARTS the first atom is the carbonyl carbon.
        acid_carbon = match[0]
        atom = mol.GetAtomWithIdx(acid_carbon)
        # Find neighbors that are carbon and not the carbonyl oxygen;
        # (in a free acid the acid carbon is bonded to one oxygen (the carbonyl) and one -OH; in an amide, it is bonded to N)
        candidate_atoms = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not candidate_atoms:
            continue
        # We expect the candidate (attachment atom) NOT to be directly in the steroid core.
        # Also, if it is not in the core, then its distance (in bonds) to the closest steroid core atom should be at least 3.
        for cand in candidate_atoms:
            if cand in steroid_core_atoms:
                # This acid group is directly fused to the ring system. Not typical for cholanic acids.
                continue
            # Calculate the shortest bond path from candidate to any atom in steroid_core_atoms.
            min_distance = 1000
            for core in steroid_core_atoms:
                sp = Chem.GetShortestPath(mol, cand, core)
                if sp and len(sp) - 1 < min_distance:
                    min_distance = len(sp) - 1
            # Require that the side chain is at least 3 bonds long (i.e. acid not directly at the core).
            if min_distance >= 3:
                sidechain_ok = True
                break
        if sidechain_ok:
            break

    if not sidechain_ok:
        return False, "Acid/conjugate group appears directly attached to the steroid core; expected an acid side chain"

    return True, "Molecule shows a fused steroid nucleus with an acid (or conjugated acid) side chain and extra hydroxyls typical for bile acids"
    
# For testing, one could call:
if __name__ == "__main__":
    # Example bile acid SMILES (one of the true positives from the reference list)
    test_smiles = "O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(C[C@H](O)CC4)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CCC(O)=O)C)[H])[H])C"
    result, reason = is_bile_acid(test_smiles)
    print(result, reason)