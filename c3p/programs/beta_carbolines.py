"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: beta-carbolines 
Definition: Any pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.
A beta‐carboline core is built on a fused bicyclic system in which a five‐membered “indole” ring 
(with an –NH group) is fused to a six‐membered ring that contains at least one nitrogen (making it a pyridine).
This method checks for a pair of fused rings that share exactly two atoms. The 5-membered ring must have exactly 
one nitrogen that carries at least one hydrogen (as expected for an indole N–H) while the 6-membered ring must contain 
at least one nitrogen.
Note: This approach does not cover every possible edge case, but is an attempt to improve over previous methods.
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string)
    appears to contain a beta-carboline scaffold or its hydrogenated derivative.
    
    Strategy:
      1. Parse the SMILES and add explicit hydrogens (to check for an –NH group).
      2. Obtain ring information.
      3. Look for a fused ring pair that share exactly two atoms where:
           - One ring is five-membered and has exactly one nitrogen that has at least one hydrogen.
           - The other ring is six-membered and contains at least one nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple (True, reason) if the scaffold is found; otherwise (False, reason).
    """
    # Parse SMILES and add explicit hydrogens so we can check N–H
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # Get rings (list of tuples of atom indices in each ring)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule."
    
    # Iterate over all pairs of rings to look for a fused system.
    for i, ring1 in enumerate(rings):
        set1 = set(ring1)
        for ring2 in rings[i+1:]:
            set2 = set(ring2)
            shared_atoms = set1.intersection(set2)
            # For a genuine fused system in a bicyclic core we expect exactly two shared atoms.
            if len(shared_atoms) != 2:
                continue
            
            # Check if one ring is 5-membered and the other is 6-membered.
            if (len(ring1) == 5 and len(ring2) == 6) or (len(ring1) == 6 and len(ring2) == 5):
                # Identify which ring is five-membered (candidate indole part) and which is six-membered (candidate pyridine part).
                ring5 = ring1 if len(ring1) == 5 else ring2
                ring6 = ring1 if len(ring1) == 6 else ring2
                
                # Count nitrogen atoms in the 5-membered ring and check for an NH.
                nitrogen_ring5 = []
                for idx in ring5:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 7:
                        nitrogen_ring5.append(atom)
                if len(nitrogen_ring5) != 1:
                    continue
                # Check that the single N in the 5-membered ring has at least one hydrogen.
                if nitrogen_ring5[0].GetTotalNumHs() < 1:
                    continue
                
                # Count nitrogen atoms in the 6-membered ring (at least one required).
                n6_count = sum(1 for idx in ring6 if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                if n6_count < 1:
                    continue
                
                # For additional confidence, demand that the two shared atoms are carbons.
                shared_ok = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in shared_atoms)
                if not shared_ok:
                    continue

                return True, ("Molecule contains a fused 5- and 6-membered ring system with a 5-membered ring "
                              "having a single NH and a 6-membered ring with at least one nitrogen (expected for a "
                              "beta-carboline or its hydrogenated derivative).")
    
    return False, ("Molecule does not appear to contain the beta-carboline scaffold based on fused ring analysis: "
                   "no fused 5- and 6-membered ring system with the required nitrogen features was found.")

# Example usage when running the script directly.
if __name__ == "__main__":
    # Some test examples (the provided list can be tested)
    test_examples = {
        "Harmine": "COc1ccc2c(c1)[nH]c1c(C)nccc12",
        "4,21-Dehydrocorynantheine aldehyde": "COC(=O)C(=C/O)\\[C@H]1C[C@H]2c3[nH]c4ccccc4c3CC[N+]2=C[C@@H]1C=C",
        "Carbazole (likely false positive)": "c1ccc2c(c1)[nH]c1ccccc12",
    }
    
    for name, smi in test_examples.items():
        result, reason = is_beta_carbolines(smi)
        print(f"{name}: {result}\n  Reason: {reason}\n")