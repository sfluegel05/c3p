"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole
Defined as: Any monocyclic heteroarene consisting of a five‐membered aromatic ring 
that contains at least one nitrogen. (The ring may also contain O or S; other elements are disallowed.)
Additional heuristics are applied to reduce false positives from large multicyclic molecules or peptides.
Heuristics:
  1. The candidate ring must be non‐fused (not sharing ≥2 atoms with any other ring).
  2. The candidate ring’s 5 atoms must comprise at least 10% of the molecule’s heavy atoms.
  3. The candidate ring must have no formally charged atoms (to avoid e.g. imidazolium salts).
  4. If the molecule has a relatively high number of amide bonds (by absolute count or ratio vs heavy atoms)
     it is assumed to be a peptide and the candidate is rejected.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_azole(smiles: str):
    """
    Determines if a molecule is (essentially) an azole based on its SMILES string.
    
    The molecule is classified as an azole if it contains a non‐fused five‐membered aromatic heterocycle 
    (with atoms limited to C, N, O or S) that contains at least one nitrogen and is a significant part 
    of the molecule (≥10% of heavy atoms). Additionally, if the molecule appears peptide‐like 
    (either ≥4 amide bonds or amide-bond:heavy-atom ratio ≥ 0.15) it is rejected.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an azole; False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"

    # Count heavy atoms
    num_heavy = mol.GetNumHeavyAtoms()
    if num_heavy == 0:
        return False, "No heavy atoms in molecule"

    # Heuristic: if there are many amide bonds relative to heavy atoms, suspect a peptide.
    # First count using SMARTS for an amide: [#6X3](=O)[#7]
    amide_smarts = Chem.MolFromSmarts("[#6X3](=O)[#7]")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    num_amide = len(amide_matches)
    # If 4 or more amide bonds OR if the ratio is high (>0.15), assume peptide.
    if num_amide >= 4 or (num_amide / num_heavy) >= 0.15:
        return False, f"Molecule contains {num_amide} amide bonds (ratio {num_amide/num_heavy:.2f}), likely a peptide"

    # Get ring information
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found in the molecule"

    # Allowed atomic numbers (C, N, O, S)
    allowed_atoms = {6, 7, 8, 16}
    candidate_rings = []
    
    # Loop over rings looking for five-membered aromatic rings that:
    #  - consist only of allowed atoms
    #  - contain at least one nitrogen
    #  - are aromatic (all atoms flagged as aromatic)
    for ring in all_rings:
        if len(ring) != 5:
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Only consider if every atom is aromatic and is allowed
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        if not all(atom.GetAtomicNum() in allowed_atoms for atom in ring_atoms):
            continue
        if not any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            continue
        
        # Reject candidate if any atom in this ring carries a non-zero formal charge.
        if any(atom.GetFormalCharge() != 0 for atom in ring_atoms):
            continue

        candidate_rings.append(ring)

    if not candidate_rings:
        return False, "No five‐membered aromatic heterocycle with nitrogen found"
    
    # Filter out fused rings: reject a candidate if it shares 2 or more atoms with any other ring.
    non_fused_candidates = []
    for cand in candidate_rings:
        cand_set = set(cand)
        is_fused = False
        for other in all_rings:
            if set(other) == cand_set:
                continue
            if len(cand_set.intersection(other)) >= 2:
                is_fused = True
                break
        if not is_fused:
            non_fused_candidates.append(cand)
    
    if not non_fused_candidates:
        return False, "No non‐fused five‐membered aromatic heterocycle with nitrogen found"
    
    # Check that at least one candidate ring is a significant portion of the molecule:
    # Since every candidate contributes 5 atoms, check if 5/num_heavy >= 0.10.
    fraction = 5 / num_heavy
    if fraction < 0.10:
        return False, f"Candidate azole ring accounts for only {fraction*100:.1f}% of heavy atoms"
    
    return True, "Found a non‐fused five‐membered aromatic heterocycle with nitrogen (azole ring)"

# Example usage (for testing; remove or comment out in production)
if __name__ == "__main__":
    test_examples = [
        # True positives
        ("CC1CCCCC1N2C(=NNC2=S)C3=CN=CC=C3", "4-(2-methylcyclohexyl)-3-(3-pyridinyl)-1H-1,2,4-triazole-5-thione"),
        ("Cn1cnc(CCN)c1", "N(tele)-methylhistamine"),
        ("CC(C)C1=C(C(=O)NC2=CC=CC=C2)C(=C(N1CC[C@@H](O)C[C@@H](O)CC(O)=O)C1=CC=C(F)C=C1)C1=CC=CC=C1", "atorvastatin"),
        ("O1C(=NC=C1)CCCCC", "2-Pentyloxazole"),
        # False positives that our filters are meant to cut out:
        ("OC(=O)[C@@H](NC(=O)CNC(=O)[C@@H](N)CC1NC=NC1)CC1NC=NC1", "Ala-His-Ser (peptide-like)"),
        ("CCCCCCCn1cc[n+](C)c1", "3-methyl-1-octylimidazolium (charged)")
    ]
    for smi, name in test_examples:
        result, explanation = is_azole(smi)
        print(f"Test: {name}\n  SMILES: {smi}\n  Result: {result}\n  Explanation: {explanation}\n")