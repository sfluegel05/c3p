"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole
Defined as: Any monocyclic heteroarene consisting of a five‐membered aromatic ring 
that contains at least one nitrogen. (The ring may also contain O or S as the only other heteroatoms.)
Heuristics added:
  1. Look for any five‐membered aromatic ring with allowed atoms (C, N, O, S), at least one N, and no formal charge.
  2. Reject if the molecule appears peptide‐like – here if it contains 2 or more amide bonds and at least 20 heavy atoms.
  3. For large molecules (≥31 heavy atoms), require that the candidate five‐membered ring represents at least 10% of the heavy atoms.
  
Note: We no longer insist that the candidate ring be “non‐fused”, so that fused rings that embed an azole core (e.g. in pyrazolopyrimidines)
are captured.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole (or contains the essential five‐membered azole ring) based on its SMILES string.
    
    The molecule is classified as an azole if it contains at least one five‐membered aromatic ring built only from
    C, N, O, and S (with no formal charge on any ring atom) that contains at least one nitrogen.
    In addition, molecules that appear peptide‐like (having ≥2 amide bonds and >20 heavy atoms) are rejected.
    For larger molecules (more than 30 heavy atoms), the azole ring is required to account for at least 10% of heavy atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as an azole; False otherwise.
        str: Explanation of the classification decision.
    """
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
    
    # Peptide check – if the molecule contains 2 or more amide bonds and is not extremely small, reject.
    # (SMARTS pattern for an amide bond: carbon with degree3 (typically in C(=O)) attached to a nitrogen)
    amide_smarts = Chem.MolFromSmarts("[#6X3](=O)[#7]")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    num_amide = len(amide_matches)
    if num_amide >= 2 and num_heavy >= 20:
        return False, f"Molecule contains {num_amide} amide bonds (with {num_heavy} heavy atoms), suggesting a peptide-like structure"
    
    # Allowed atomic numbers: Carbon (6), Nitrogen (7), Oxygen (8), Sulfur (16)
    allowed_atoms = {6, 7, 8, 16}
    
    # Look for candidate rings among all rings reported by RDKit.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    candidate_found = False
    candidate_details = ""
    
    # Loop over every ring; we consider five-membered rings only.
    for ring in atom_rings:
        if len(ring) != 5:
            continue
        # Get atoms in this ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check that all atoms are aromatic
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        # Check that every atom is one of the allowed types (C, N, O, S)
        if not all(atom.GetAtomicNum() in allowed_atoms for atom in ring_atoms):
            continue
        # There must be at least one nitrogen in the ring
        if not any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            continue
        # Reject if any ring atom carries a nonzero formal charge
        if any(atom.GetFormalCharge() != 0 for atom in ring_atoms):
            continue
        # If we get here, this ring is a candidate azole ring.
        candidate_found = True
        # (We break out after the first valid candidate is found.)
        break
    
    if not candidate_found:
        return False, "No five‐membered aromatic heterocycle with at least one nitrogen (azole ring) found"
    
    # For small molecules the candidate ring is more significant.
    # For larger molecules, require that five atoms represent at least 10% of heavy atoms.
    if num_heavy > 30:
        fraction = 5 / num_heavy
        if fraction < 0.10:
            return False, f"Candidate azole ring accounts for only {fraction*100:.1f}% of heavy atoms (molecule too large)"
    
    return True, "Found a five‐membered aromatic heterocycle with at least one nitrogen (azole ring) that meets the criteria"

# Example usage (for testing; remove or comment out in production)
if __name__ == "__main__":
    test_examples = [
        # True positives
        ("CC1CCCCC1N2C(=NNC2=S)C3=CN=CC=C3", "4-(2-methylcyclohexyl)-3-(3-pyridinyl)-1H-1,2,4-triazole-5-thione"),
        ("Cn1cnc(CCN)c1", "N(tele)-methylhistamine"),
        ("CC(C)C1=C(C(=O)NC2=CC=CC=C2)C(=C(N1CC[C@@H](O)C[C@@H](O)CC(O)=O)C1=CC=C(F)C=C1)C1=CC=CC=C1", "atorvastatin"),
        ("CC(C)(C)C(O)C(=C/C1CCCCC1)\n1cncn1", "Triapenthenol"),
        ("COC1=CC=CC=C1C=NNC(=O)C2=C(N(N=N2)C3=NON=C3N)C4=CC=CC=C4", "1-(4-amino-1,2,5-oxadiazol-3-yl)-N-[(2-methoxyphenyl)methylideneamino]-5-phenyl-4-triazolecarboxamide"),
        # False positives (expected rejection)
        ("OC(=O)[C@@H](NC(=O)CNC(=O)[C@@H](N)CC1NC=NC1)CC1NC=NC1", "Ala-His-Ser (peptide-like)"),
        ("CN1C=C(N=C1)S(=O)(=O)N2CC[C@H]3[C@@H]2C4=C(C=CC(=C4)C5=CC(=CC=C5)OC)N[C@H]3CO", "Charged/fused system")
    ]
    
    for smi, name in test_examples:
        result, explanation = is_azole(smi)
        print(f"Test: {name}\n  SMILES: {smi}\n  Result: {result}\n  Explanation: {explanation}\n")