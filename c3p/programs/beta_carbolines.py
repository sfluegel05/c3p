"""
Classifies: CHEBI:60834 beta-carbolines
"""
#!/usr/bin/env python
"""
Classifies: beta-carbolines 
Definition: Any pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.
A beta‐carboline core is built on a fused bicyclic system in which a five‐membered “indole” ring 
(with an –NH group) is fused to a six‐membered ring that contains at least one nitrogen (making it a pyridine).
This approach analyzes the rings present in a molecule and checks for a fused 5– and 6–membered system with these features.
Note: This method is an approximation and may still mis‐classify some edge cases.
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) 
    appears to contain a beta-carboline (or hydrogenated beta-carboline) scaffold.
    
    The strategy is:
      1. Parse the SMILES into an RDKit molecule.
      2. Obtain the ring information.
      3. Look for a fused ring pair where one ring is five-membered containing exactly one nitrogen
         (the “indole” part) and the other is six-membered containing at least one nitrogen (the “pyridine” part).
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple (True, reason) if the scaffold is found; otherwise (False, reason).
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get ring information (list of tuples; each tuple contains indices of atoms in the ring)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule."
    
    # Iterate through all pairs of rings and check for fused rings sharing at least 2 atoms.
    for i, ring1 in enumerate(rings):
        set1 = set(ring1)
        for ring2 in rings[i+1:]:
            set2 = set(ring2)
            shared_atoms = set1.intersection(set2)
            if len(shared_atoms) < 2:
                # not a fused ring system
                continue
            
            # Look for a pair where one ring is 5-membered and the other is 6-membered.
            if (len(ring1) == 5 and len(ring2) == 6) or (len(ring1) == 6 and len(ring2) == 5):
                # Identify which ring is 5-membered (expected to be the indole part) and which is 6-membered (expected pyridine)
                ring5 = ring1 if len(ring1) == 5 else ring2
                ring6 = ring1 if len(ring1) == 6 else ring2
                
                # Count nitrogen atoms in the 5-membered ring.
                n_count_ring5 = sum(1 for idx in ring5 if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                # For a proper indole-like ring, we expect exactly one nitrogen.
                if n_count_ring5 != 1:
                    continue
                
                # Count nitrogen atoms in the 6-membered ring.
                n_count_ring6 = sum(1 for idx in ring6 if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                if n_count_ring6 < 1:
                    # This excludes fused systems with a benzene ring instead of a pyridine.
                    continue

                # Additional check: you might require that the nitrogen in the five-membered ring has at least one hydrogen,
                # but in hydrogenated variants the explicit hydrogen may not be present. So here we simply require one nitrogen.
                
                return True, ("Molecule contains a fused 5- and 6-membered ring system with appropriate nitrogen counts "
                              "(expected for a beta-carboline or its hydrogenated derivative).")
    
    return False, "Molecule does not appear to contain the beta-carboline scaffold (a fused indole-pyridine system) based on ring analysis."

# Example usage (for testing) when running the script directly.
if __name__ == "__main__":
    test_examples = {
        "Aromatic harmine": "COc1ccc2c(c1)[nH]c1cnccc12",  # harmine core (modified to force pyridine N)
        "Hydrogenated variant": "O(C1=CC=2C3=C(NC2C=C1)CN(CC3)C)",  # one of the provided examples
        "Carbazole (should be false positive)": "c1ccc2c(c1)[nH]c1ccccc12",
    }
    
    for name, smi in test_examples.items():
        result, reason = is_beta_carbolines(smi)
        print(f"{name}: {result}\n  Reason: {reason}\n")