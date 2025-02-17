"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent.
A valid molecule must have exactly one carboxylic acid group 
(–C(=O)[OH] or –C(=O)[O-]) and exactly one α‐carbon (the carbon directly 
bonded to the carboxyl carbon) that bears a 2‑oxo substituent.
This implementation checks if the α‐carbon (or a substituent on the α‐carbon) 
has a C=O bond.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    It first verifies that there is exactly one carboxylic acid group.
    Then, for the unique acid group, it checks its directly bonded carbon neighbors
    (α‐carbons) for evidence of a 2‑oxo substituent. The oxo substituent may be present in two forms:
      Option A: the α‐carbon itself is a carbonyl (has C(=O) as part of its bonds);
      Option B: the α‐carbon is sp3 but carries a substituent carbon that is carbonyl (has a double bond to oxygen).
      
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): Tuple with a boolean (True if classification criteria met) 
                     and a string explaining the reason.
                     In case of a parsing error, returns (False, "Invalid SMILES string").
    """
    
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for a carboxylic acid group (covers –COOH and –COO– forms)
    acid_smarts = "[CX3](=O)[OX2H1,O-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Must have exactly one carboxylic acid group.
    if len(acid_matches) == 0:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups; expected exactly one (monocarboxylic acid)"
    
    # For the unique hit, assume the first atom in the SMARTS (index 0) is the acid carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Get all neighboring atoms of the acid carbon that are carbons (potential α‐carbons).
    alpha_candidates = []
    for nbr in acid_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # must be carbon
            alpha_candidates.append(nbr)
    
    if not alpha_candidates:
        return False, ("No alpha carbon found. Likely a formic acid derivative so there is no "
                       "possibility for a 2‑oxo substituent")
                       
    valid_alpha = []
    
    # Function to check if a given atom has a double bond to oxygen.
    def has_carbonyl(atom):
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:  # oxygen
                    return True
        return False

    # Evaluate each alpha candidate using two possible options:
    # Option A: The candidate itself displays a carbonyl (C=O).
    # Option B: The candidate (which might be sp3) has a substituent (other than the acid carbon)
    #           that behaves as a carbonyl (i.e. has a C=O).
    for alpha in alpha_candidates:
        qualifies = False
        # Option A
        if has_carbonyl(alpha):
            qualifies = True
        else:
            # Option B: check substituents of alpha (excluding the acid carbon)
            for sub in alpha.GetNeighbors():
                if sub.GetIdx() == acid_carbon_idx:
                    continue
                if sub.GetAtomicNum() == 6 and has_carbonyl(sub):
                    qualifies = True
                    break
        if qualifies:
            valid_alpha.append(alpha)
    
    if len(valid_alpha) == 0:
        return False, ("No α‐carbon with an attached 2‑oxo (C=O) substituent found; "
                       "molecule lacks a clear 2‑oxo group")
    if len(valid_alpha) > 1:
        return False, (f"Found {len(valid_alpha)} potential α‐carbons with a 2‑oxo substituent — "
                       "ambiguous for classification")
    
    return True, ("Contains a single carboxylic acid group with an α‐carbon (or attached substituent) "
                  "bearing a 2‑oxo (C=O) group")

# Example usage:
# Uncomment for testing:
# test_smiles = "CCCCCCCC(=O)C(O)=O"  # 2-oxononanoic acid (a true positive)
# result, reason = is_2_oxo_monocarboxylic_acid(test_smiles)
# print(result, reason)