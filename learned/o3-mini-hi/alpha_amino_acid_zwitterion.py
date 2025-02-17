"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
#!/usr/bin/env python
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid–zwitterion obtained by transfer of a proton from the carboxyl 
group to the amino group of any alpha–amino acid; major species at pH 7.3.
This implementation first finds a carboxylate group (C(=O)[O-]) and then checks 
if its attached (alpha) carbon is bonded to a protonated (positively charged) nitrogen.
"""

from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines whether a molecule represents an alpha-amino-acid zwitterion.
    The heuristic is as follows:
       1. Identify the carboxylate group using the SMARTS "C(=O)[O-]".
       2. For each carboxylate group found, check its connected atom (the candidate alpha carbon).
       3. At the candidate alpha carbon, verify whether it has at least one neighboring nitrogen
          with a positive formal charge (e.g. [NH3+] or [NH+]). This corresponds to the
          protonated amino group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if an alpha-amino-acid zwitterion pattern is found, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to help with valence/hydrogen count if needed.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS for a carboxylate group: C(=O)[O-]
    carboxylate_smarts = "C(=O)[O-]"
    carboxylate_query = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_query is None:
        return False, "Error in carboxylate SMARTS pattern definition"
    
    # Find all matches for the carboxylate group
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_query)
    if not carboxylate_matches:
        return False, "No carboxylate (C(=O)[O-]) group found"
    
    # For each carboxylate match, check if the carboxyl carbon is attached to an alpha carbon 
    # that itself is bonded to a protonated nitrogen.
    for match in carboxylate_matches:
        # In the matched pattern, we expect:
        # match[0]: carboxyl carbon, match[1]: the carbonyl oxygen, match[2]: the deprotonated oxygen.
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
        # Iterate over neighbors of the carboxyl carbon to find candidate alpha carbon(s)
        for neighbor in carboxyl_carbon.GetNeighbors():
            # Skip oxygen atoms from the carboxylate group (we want the carbon chain)
            if neighbor.GetAtomicNum() == 8:
                continue
            # The candidate alpha carbon must be a carbon (typically sp3, possibly chiral)
            if neighbor.GetAtomicNum() != 6:
                continue
            alpha_carbon = neighbor
            # Now, check if alpha carbon has any neighbor that is a protonated nitrogen.
            # (A protonated nitrogen carries a positive charge; typical representations are [NH3+]
            # or [NH+]. We check if its formal charge is positive.)
            for nb in alpha_carbon.GetNeighbors():
                if nb.GetAtomicNum() == 7 and nb.GetFormalCharge() > 0:
                    # We found a candidate alpha carbon connected to a positive nitrogen.
                    # For the reason string we report the indices.
                    return True, ("Alpha-amino acid zwitterion pattern found: "
                                  "carboxylate group on atom idx {} is attached to alpha carbon (idx {}) "
                                  "which is bonded to a protonated nitrogen (idx {})."
                                  .format(carboxyl_carbon_idx, alpha_carbon.GetIdx(), nb.GetIdx()))
    return False, ("No alpha-amino acid zwitterion backbone detected: "
                   "could not find a carboxylate group whose attached carbon is also bonded "
                   "to a protonated nitrogen.")

# Example usage with some test SMILES from the provided list.
if __name__ == "__main__":
    test_smiles = [
        "[C@H]1(CCC[NH+]1C)C([O-])=O",                 # N-methylproline zwitterion
        "[C@@H]([C@@H](C([O-])=O)[NH3+])(O)C#C",         # L-beta-ethynylserine zwitterion
        "[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)[O-])=O)NC([C@@H](C3=CC=C(C=C3)OCC[C@H](C([O-])=O)[NH3+])[NH3+])=O",  # nocardicin C dizwitterion
        "CC[C@H](C)[C@H]([NH3+])C([O-])=O",              # L-isoleucine zwitterion
        "CCCC[C@H]([NH3+])C([O-])=O",                    # L-2-aminohexanoic acid zwitterion
        "[NH3+][C@@H](CO)C([O-])=O",                     # L-serine zwitterion
    ]
    
    for smi in test_smiles:
        result, reason = is_alpha_amino_acid_zwitterion(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)