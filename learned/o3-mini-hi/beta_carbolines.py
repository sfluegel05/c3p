"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: beta-carbolines
Definition: Any pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.
A beta‐carboline core is built on a fused bicyclic system where a five‐membered “indole” ring 
(with an –NH group) is fused to a six‐membered ring that contains at least one nitrogen (making it a pyridine).
This method attempts to identify such a pair by:
  1. Finding all 5- and 6-membered rings (using GetSymmSSSR).
  2. Requiring that a candidate 5-membered ring contains exactly one nitrogen that bears at least one hydrogen.
  3. Requiring that a candidate 6-membered ring contains at least one nitrogen.
  4. Looking for a pair of such rings that share exactly two atoms (consistent with a fused bicyclic system).
Note: This heuristic does not capture every edge case but is an improved approach over the previous code.
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string)
    appears to contain a beta-carboline scaffold (or its hydrogenated derivative).

    Strategy:
      1. Parse the SMILES and add explicit hydrogens (to check for an –NH group).
      2. Compute the smallest set of smallest rings (SSSR) and filter for 5- and 6-membered rings.
      3. Require that one of the 5-membered rings (candidate indole) has exactly one nitrogen with ≥1 hydrogen.
      4. Require that one of the 6-membered rings (candidate pyridine) has at least one nitrogen.
      5. Check if the two rings are fused, meaning they share exactly two atoms.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): Tuple; True with a reason if the scaffold is found, otherwise False with a reason.
    """
    # Parse the molecule and add hydrogens so we can check for N–H.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # Obtain the smallest set of smallest rings.
    rings = list(Chem.GetSymmSSSR(mol))
    if not rings:
        return False, "No rings found in the molecule."
    
    # Filter for 5-membered and 6-membered rings.
    rings5 = [set(r) for r in rings if len(r) == 5]
    rings6 = [set(r) for r in rings if len(r) == 6]
    if not rings5:
        return False, "No 5-membered rings found."
    if not rings6:
        return False, "No 6-membered rings found."
    
    # Iterate over candidate 5-membered rings (indole part)
    for r5 in rings5:
        # Count nitrogen atoms in the candidate 5-membered ring.
        nitrogen_indices = [idx for idx in r5 if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7]
        if len(nitrogen_indices) != 1:
            continue  # must have exactly one nitrogen
        # Ensure that the nitrogen has at least one hydrogen (as in an NH group).
        n_atom = mol.GetAtomWithIdx(nitrogen_indices[0])
        if n_atom.GetTotalNumHs() < 1:
            continue
        
        # Now iterate over 6-membered ring candidates (pyridine part).
        for r6 in rings6:
            # Check if the two rings share exactly 2 atoms (typical for a fused bicyclic system).
            shared_atoms = r5.intersection(r6)
            if len(shared_atoms) != 2:
                continue
            
            # Check that the 6-membered ring has at least one nitrogen.
            if sum(1 for idx in r6 if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7) < 1:
                continue
            
            # If all conditions are met, we assume the beta-carboline scaffold is present.
            return True, ("Molecule contains a fused 5-membered ring (indole with an NH) and a fused "
                          "6-membered ring with a nitrogen, consistent with a beta-carboline scaffold or its derivative.")
    
    # If no valid fused ring pair is found.
    return False, ("Molecule does not appear to contain the beta-carboline scaffold based on fused ring analysis: "
                   "no fused 5- and 6-membered ring system with the required nitrogen features was found.")

# Example usage: run some tests when executing the script directly.
if __name__ == "__main__":
    # Test some example SMILES (some from the provided list)
    test_examples = {
        "Harmine": "COc1ccc2c(c1)[nH]c1c(C)nccc12",
        "4,21-Dehydrocorynantheine aldehyde": "COC(=O)C(=C/O)\\[C@H]1C[C@H]2c3[nH]c4ccccc4c3CC[N+]2=C[C@@H]1C=C",
        "Carbazole (likely negative)": "c1ccc2c(c1)[nH]c1ccccc12",
    }
    for name, smi in test_examples.items():
        result, reason = is_beta_carbolines(smi)
        print(f"{name}: {result}\n  Reason: {reason}\n")