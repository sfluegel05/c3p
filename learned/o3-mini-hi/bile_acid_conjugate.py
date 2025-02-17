"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile acid conjugates.
A bile acid conjugate is defined as any bile acid (a molecule with a steroid nucleus --
typically a fused 4‐ring structure with at least one five–membered ring and one six–membered ring,
and a carboxylate side chain that has been converted to an amide or ester) that is conjugated
to a functional group that increases hydrophilicity or charge (for example glycine, taurine,
sulfate, glucuronate, or an uncharged sugar).
This implementation uses heuristic SMARTS‐based filters using rdkit.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.

    This function uses two main heuristics:
    1. The molecule should contain a “steroid‐like” bile acid core. We require that the
       molecule has at least 4 rings and that among these rings at least one is five–membered
       and one is six–membered. (This is a rough filter for the cholane nucleus.)
    2. The molecule should show evidence of conjugation. In conjugated bile acids the carboxyl 
       is converted to an amide or ester linkage. (a) We require the presence of at least one 
       such acyclic amide or ester bond; and (b) we search for one of several conjugate fragment
       patterns (e.g. glycine, taurine, sulfate, glucuronate or sugar).
       
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a bile acid conjugate, False otherwise.
       str: Reason explaining the classification decision.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ----- Heuristic 1: Look for steroid (bile acid) nucleus ----- #
    # Instead of just count rings, we check:
    #   (a) that the molecule contains at least 4 rings overall,
    #   (b) that at least one ring is five-membered and at least one is six-membered.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) < 4:
        return False, f"Insufficient rings: found {len(atom_rings)}, expected at least 4 for a bile acid core."
    
    ring_sizes = [len(ring) for ring in atom_rings]
    if 5 not in ring_sizes:
        return False, "No five-membered ring found (expected one D-ring in a bile acid nucleus)."
    if 6 not in ring_sizes:
        return False, "No six-membered ring found (expected at least one six-membered ring in the steroid nucleus)."
    
    # ----- Heuristic 2a: Confirm presence of an amide or ester bond (non-ring) ----- #
    # A conjugated bile acid should have a linkage from its side chain (converted carboxyl group).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ester_pattern   = Chem.MolFromSmarts("C(=O)O")
    linkage_found = False
    for pattern in (amide_pattern, ester_pattern):
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        # Check if at least one match is on an acyclic (non-ring) bond.
        for match in matches:
            # match is a tuple of atom indices corresponding to C, =O, and N or O
            # Here we check that the carbonyl carbon is not in a ring:
            if not mol.GetAtomWithIdx(match[0]).IsInRing():
                linkage_found = True
                break
        if linkage_found:
            break
    if not linkage_found:
        return False, "No acyclic amide/ester linkage found that is expected for a conjugated bile acid."

    # ----- Heuristic 2b: Look for one of several conjugate fragments ----- #
    # We set up a list of conjugate SMARTS patterns.
    # (Note: these patterns are simplified and may not cover all edge cases.)
    conjugate_patterns = [
        # glycine-like: look for an amide-bound CH2-COOH (optionally charged)
        ("glycine", "[NX3;H2,H1][CH2]C(=O)[O-]?"),
        # taurine: a 2-aminoethanesulfonic acid; note the pattern may allow optional protonation
        ("taurine", "NCCS(=O)(=O)[O-]?"),
        # sulfate: an OSO3– pattern
        ("sulfate", "[OS](=O)(=O)[O-]?"),
        # glucuronate: simplified pattern for glucuronic acid fragment attached via an O-linkage.
        ("glucuronate", "OC1C(=O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C1"),
        # sugar: simplified pattern for a hexose (may be noisy; usually in conjugation one sees an ester/amide linkage too)
        ("sugar", "OC1C(O)C(O)C(O)C(O)C1O")
    ]
    
    conjugate_found = False
    conjugate_reason = ""
    for name, smarts in conjugate_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # skip if SMARTS is invalid
        if mol.HasSubstructMatch(pattern):
            conjugate_found = True
            conjugate_reason = f"Found conjugate group pattern matching {name}."
            break

    if not conjugate_found:
        return False, "No conjugate fragment detected (expected glycine, taurine, sulfate, glucuronate, or sugar)."
    
    # ----- Final decision ----- #
    # If both the steroid nucleus (with characteristic ring sizes) and the linkage plus conjugate fragment are present,
    # then classify as a bile acid conjugate.
    return True, f"Contains a steroid nucleus with sufficient rings and a conjugate fragment. {conjugate_reason}"

# Example usage (for testing):
if __name__ == "__main__":
    # A known bile acid conjugate (one of the provided examples)
    test_smiles = "O[C@H]1[C@]2([C@]([C@]3([C@@]([C@](CC3)([C@@H](CCC(=O)NC(CO)C(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(C[C@H](O)CC4)[H])C)[H])C)[H])[H]"
    result, reason = is_bile_acid_conjugate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)