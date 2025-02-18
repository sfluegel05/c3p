"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Hydroxy-5β-cholanic acids (bile acids)
Definition: Bile acids are hydroxy-5β-cholanic acids (or their glycine/taurine conjugates). They typically
display a fused steroid nucleus (one five-membered ring fused with three six-membered rings), an acid group (or an amide
formed with glycine/taurine) that is located on an acid side chain (i.e. several bonds away from the poly‐ring core),
and at least one extra hydroxyl substituent.
The following implementation uses several SMARTS and distance/regional checks to improve specificity.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid (a hydroxy-5β-cholanic acid or its glycine/taurine conjugate)
    based on its SMILES string.
    
    The molecule must pass these tests:
      1. Be parsed correctly.
      2. Contain either a free carboxylic acid group OR show an amide bond arising from conjugation with a glycine or taurine fragment.
         (Three SMARTS patterns are used.)
      3. Possess at least three hydroxyl groups ([OX2H]) (one may be part of the acid motif but an extra OH is expected).
      4. Have a molecular weight in roughly the 350–700 Da range.
      5. Contain a fused steroid nucleus: we require at least one five-membered ring and at least three six-membered rings.
         Moreover, we try to “assemble” a fused core by gathering atoms from rings of size 5 or 6 and require a minimum size.
      6. The acid (or conjugate) group should appear on a side chain. Specifically, the carbonyl carbon (of the acid motif)
         must be on a non‐ring atom whose attachment to the steroid core (determined by a bond shortest-path) is at least 4 bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a bile acid, False otherwise.
        str: A reason giving the basis for the classification decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the acid/conjugate group.
    # Free carboxylic acid: carbonyl carbon not in a ring, with an –OH.
    acid_pattern = Chem.MolFromSmarts("[$([CX3](=O)[O;H,-]);!R]")
    # Glycine conjugate: carbonyl carbon with an amide linked to –CH2C(=O)[O;H,-]
    glycine_pattern = Chem.MolFromSmarts("[$([CX3](=O)NCC(=O)[O;H,-]);!R]")
    # Taurine conjugate: carbonyl carbon with an amide linked to –CH2CH2S(=O)(=O)[O;H,-]
    taurine_pattern = Chem.MolFromSmarts("[$([CX3](=O)NCCS(=O)(=O)[O;H,-]);!R]")

    acid_matches = mol.GetSubstructMatches(acid_pattern)
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    taurine_matches = mol.GetSubstructMatches(taurine_pattern)
    
    if not (acid_matches or glycine_matches or taurine_matches):
        return False, "Missing free acid or recognized conjugated acid group essential for cholanic acids"
    
    # Check for hydroxyl groups ([OX2H]) and require at least 3 total.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Too few hydroxyl substituents (found {len(hydroxyl_matches)}; expect at least 3)"
    
    # Check molecular weight in the 350-700 Da range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 700):
        return False, f"Molecular weight {mol_wt:.1f} out of expected bile acid range (350-700 Da)"
    
    # Analyze ring systems and fused steroid nucleus.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found; expected a steroid nucleus"
    
    # Count number of 5- and 6-membered rings (these are typical in a steroid nucleus).
    count_5 = sum(1 for ring in rings if len(ring) == 5)
    count_6 = sum(1 for ring in rings if len(ring) == 6)
    if count_5 < 1 or count_6 < 3:
        return False, f"Ring system does not meet steroid criteria (found {count_5} 5-membered and {count_6} 6-membered rings)"
    
    # Collate all atoms that appear in any 5- or 6-membered ring. We expect the fused core to be large (at least 15 atoms).
    steroid_core_atoms = set()
    for ring in rings:
        if len(ring) in (5,6):
            for idx in ring:
                steroid_core_atoms.add(idx)
    if len(steroid_core_atoms) < 15:
        return False, "Fused steroid core not clearly discernible (not enough atoms in 5/6 membered rings)"
    
    # Next, check that the acid group is on a side-chain away from the steroid core.
    # We combine matches from free acid and conjugate SMARTS and require that the attachment (neighbor carbon of C=O) is not in the core,
    # plus its distance to any atom in the steroid core is at least 4 bonds.
    combined_matches = list(acid_matches) + list(glycine_matches) + list(taurine_matches)
    sidechain_ok = False
    for match in combined_matches:
        # The SMARTS pattern is designed so that the first match atom is the carbonyl carbon.
        acid_carbon = match[0]
        atom_obj = mol.GetAtomWithIdx(acid_carbon)
        # Ensure the acid carbon is not itself in the steroid core.
        if acid_carbon in steroid_core_atoms:
            continue
        # Get neighboring carbon atoms (exclude heteroatoms like oxygen or N) that might be the attachment.
        cand_neighbors = [nbr.GetIdx() for nbr in atom_obj.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not cand_neighbors:
            continue
        for cand in cand_neighbors:
            # We insist the candidate attachment atom is not in the steroid core.
            if cand in steroid_core_atoms:
                continue
            # Compute the shortest path from this candidate to any atom of the steroid core.
            min_distance = 1000
            for core_atom in steroid_core_atoms:
                sp = Chem.GetShortestPath(mol, cand, core_atom)
                if sp and (len(sp) - 1) < min_distance:
                    min_distance = len(sp) - 1
            # Require the side chain to be at least 4 bonds long.
            if min_distance >= 4:
                sidechain_ok = True
                break
        if sidechain_ok:
            break

    if not sidechain_ok:
        return False, "Acid/conjugate group appears too close to the steroid core; expected a side chain at least 4 bonds long"
    
    return True, "Molecule shows a fused steroid nucleus with a separated acid (or conjugated amide) side chain and extra hydroxyls typical for bile acids"

# For testing:
if __name__ == "__main__":
    # Example: one of the true positive bile acid SMILES from the provided list.
    test_smiles = "O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(C[C@H](O)CC4)[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CCC(O)=O)C)[H])[H])C"
    result, reason = is_bile_acid(test_smiles)
    print(result, reason)