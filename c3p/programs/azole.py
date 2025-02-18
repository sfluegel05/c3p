"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole
Defined as: any monocyclic heteroarene consisting of a five‐membered aromatic ring 
that contains at least one nitrogen. (The ring may include O or S but no other elements.)
To reduce false positives from large multicyclic molecules or peptides,
this version requires that the candidate five‐membered ring is “non‐fused”
(i.e. not sharing ≥2 atoms with another ring) and that it makes up at least 10% of the molecule’s heavy atoms.
Additionally, if the molecule contains ≥4 amide bonds it is assumed to be a peptide,
and the candidate azole is rejected.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_azole(smiles: str):
    """
    Determines if a molecule is (essentially) an azole based on its SMILES string.
    
    Our algorithm checks for a five‐membered aromatic heterocycle (with C, N, O, S only),
    containing at least one nitrogen. In order to filter out false positives that are large molecules
    (in which an azole fragment would be only a small substituent) or peptides (with many amide bonds),
    we require that:
      1. The candidate ring is non‐fused (not sharing ≥2 atoms with any other ring)
      2. The five‐membered ring accounts for at least 10% of the heavy atoms.
      3. The molecule does not contain 4 or more amide bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an azole, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"

    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Allowed atomic numbers: C=6, N=7, O=8, S=16
    allowed_atoms = {6, 7, 8, 16}
    candidate_rings = []
    # First, find rings that are exactly five atoms, aromatic,
    # contain only allowed elements, and have at least one nitrogen.
    for ring in all_rings:
        if len(ring) != 5:
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        if not all(atom.GetAtomicNum() in allowed_atoms for atom in ring_atoms):
            continue
        if not any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            continue
        candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No five‐membered aromatic heterocycle with nitrogen found"
    
    # Now discard rings that are fused – i.e. that share at least 2 atoms with any other ring.
    non_fused_candidates = []
    for cand in candidate_rings:
        cand_set = set(cand)
        fused = False
        for other in all_rings:
            # Skip if same ring
            if set(other) == cand_set:
                continue
            # If the candidate shares 2 or more atoms with any other ring, we treat it as fused.
            if len(cand_set.intersection(other)) >= 2:
                fused = True
                break
        if not fused:
            non_fused_candidates.append(cand)
    
    if not non_fused_candidates:
        return False, "No non‐fused five‐membered aromatic heterocycle with nitrogen found"
    
    # Check that the candidate ring is not a very minor substructure.
    num_heavy = mol.GetNumHeavyAtoms()
    # Since every candidate here is 5 atoms, compute the fraction.
    max_fraction = 5 / num_heavy
    if max_fraction < 0.10:
        return False, f"Candidate azole ring accounts for only {max_fraction*100:.1f}% of the molecule's heavy atoms"
    
    # Heuristic to reduce peptides: count amide bonds.
    # Use a SMARTS pattern for an amide: [CX3](=O)N
    amide_smarts = Chem.MolFromSmarts("[CX3](=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 4:
        return False, f"Molecule contains {len(amide_matches)} amide bonds, likely a peptide"
    
    return True, "Found a non‐fused five‐membered aromatic heterocycle with nitrogen (azole ring)"

# Example usage (for testing; remove or comment out in production):
if __name__ == "__main__":
    examples = [
        ("CC(C)OCCCNC1=NC=NC2=C1C=NN2C3=CC=CC=C3", "1-phenyl-N-(3-propan-2-yloxypropyl)-4-pyrazolo[3,4-d]pyrimidinamine"),
        ("CC1CCCCC1N2C(=NNC2=S)C3=CN=CC=C3", "4-(2-methylcyclohexyl)-3-(3-pyridinyl)-1H-1,2,4-triazole-5-thione"),
        ("Cn1cnc(CCN)c1", "N(tele)-methylhistamine"),
        ("OC1=CC=CC=C1C=NNC(=O)C2=C(N(N=N2)C3=NON=C3N)C4=CC=CC=C4", "1-(4-amino-1,2,5-oxadiazol-3-yl)-..."),
        # A couple of examples from the false negatives:
        ("CC(C)C1=C(C(=O)NC2=CC=CC=C2)C(=C(N1CC[C@@H](O)C[C@@H](O)CC(O)=O)C1=CC=C(F)C=C1)C1=CC=CC=C1", "atorvastatin"),
        ("O=C1NC(=O)C(c2c[nH]c3ccccc23)=C1c1c[nH]cn1", "didemnimide A"),
    ]
    for smi, name in examples:
        res, reason = is_azole(smi)
        print(f"Test: {name}\n  SMILES: {smi}\n  Result: {res}\n  Reason: {reason}\n")