"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid‐zwitterion obtained by transfer of a proton from the carboxy to the amino group
of any alpha‐amino acid; major species at pH 7.3.
This improved implementation first searches for deprotonated carboxylate groups (C(=O)[O–]),
then examines each neighboring carbon (“alpha carbon”) to see whether it is sp3,
has at least one hydrogen, and is directly bonded to a positively charged nitrogen.
This two‐step process is intended to catch cases with chiral specifications ([C@H] or [C@@H])
as well as cyclic amino acids, and to help avoid some false positives.
"""

from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    
    Instead of a single SMARTS match, the method first looks for a deprotonated carboxylate group,
    then checks its neighboring carbon (the candidate alpha carbon) for being sp3 (and having at least one hydrogen)
    and for being directly bonded to a nitrogen with a positive formal charge (e.g. [NH3+] for a free amino acid
    or [NH2+] in a cyclic system such as proline derivatives).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an alpha-amino-acid zwitterion, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a deprotonated carboxylate group: C(=O)[O-]
    carboxylate_smarts = "C(=O)[O-]"
    carboxylate_pat = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_pat is None:
        return False, "Failed to create carboxylate SMARTS pattern"
    
    # Get all matches for the carboxylate pattern.
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pat)
    if not carboxylate_matches:
        return False, "No deprotonated carboxylate group (C(=O)[O-]) found"

    # For each carboxylate found, look at the carbon (first atom in our pattern) as the carboxyl carbon.
    for match in carboxylate_matches:
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
        
        # Examine each neighbor of the carboxyl carbon (ignoring the oxygens in the acid group)
        for nbr in carboxyl_carbon.GetNeighbors():
            # We assume the alpha carbon candidate is an aliphatic carbon
            if nbr.GetAtomicNum() != 6:
                continue

            # Check that the candidate alpha carbon has sp3 hybridization (allowing for chiral markers like [C@H])
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            
            # Check that the candidate alpha carbon carries at least one hydrogen
            if nbr.GetTotalNumHs() < 1:
                continue

            # Now check if this candidate is directly bonded to a nitrogen that is positively charged.
            for nbr2 in nbr.GetNeighbors():
                # Look for nitrogen (atomic number 7) with a positive formal charge.
                if nbr2.GetAtomicNum() == 7 and nbr2.GetFormalCharge() > 0:
                    # Optionally, require that the nitrogen has at least one hydrogen (to discount quaternary nitrogens).
                    if nbr2.GetTotalNumHs() < 1:
                        continue
                    # Found a candidate alpha carbon with both the -COO- and the protonated amine.
                    return True, ("Found matching alpha carbon (atom idx {}) attached to deprotonated "
                                  "carboxylate and protonated amine (atom idx {})."
                                  .format(nbr.GetIdx(), nbr2.GetIdx()))
    
    # If no alpha carbon with the full signature is found, the molecule is not classified correctly.
    return False, "No alpha carbon with both connected deprotonated carboxylate and protonated amine found"

# Example usage (you can remove or comment out the test cases):
if __name__ == "__main__":
    # Example: N-methylproline zwitterion
    test_smiles = "[C@H]1(CCC[NH+]1C)C([O-])=O"
    result, reason = is_alpha_amino_acid_zwitterion(test_smiles)
    print("SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)