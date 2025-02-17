"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: Glucosinolate
Definition: Water-soluble anionic substituted thioglucosides.
Glucosinolates have a central C atom which is bonded via an S atom to a glycone group
and via an N atom to a sulfonated oxime group, and which also carries a side‐group.
The side-chain and sulfate group have an anti stereochemical configuration across the C=N double bond.
(Note: The anti configuration is not explicitly verified.)
"""

from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.

    We require that the molecule contains:
      1. A sulfonated oxime group. We use the SMARTS:
             C(=N[OX2]S(=O)(=O)[O-]?)
         (ignoring chirality). Here the central carbon is double bonded to an N,
         which is connected to an oxygen and then to a sulfonyl group.
      2. A central carbon (from the oxime match) that has at least three substituents,
         one of which is a sulfur atom that is not the sulfonyl S from the oxime.
      3. That sulfur atom (the “glycoside S”) is further connected (by any bond) to a sugar ring.
         For the sugar (glycone) we use a relaxed SMARTS for a pyranose ring:
             C1OC(C(O)C(O)C1O)
         and we ignore chirality.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as glucosinolate, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- Step 1: Find the sulfonated oxime group ---
    # Our SMARTS pattern for the sulfonated oxime group (ignoring chirality)
    # Pattern: a carbon double bonded to an N; that N is single-bonded to an O, 
    # which is bonded to a sulfonyl group S(=O)(=O)[O-] (the last oxygen is optionally charged).
    oxime_smarts = "C(=N[OX2]S(=O)(=O)[O-]?)"
    oxime_query = Chem.MolFromSmarts(oxime_smarts)
    if oxime_query is None:
        return None, None  # unexpected error with SMARTS
    
    # Use useChirality=False to allow matches regardless of chiral annotations.
    oxime_matches = mol.GetSubstructMatches(oxime_query, useChirality=False)
    if not oxime_matches:
        return False, "No sulfonated oxime substructure found."
    
    # --- Step 2: Find sugar (glycone) ring substructure ---
    # We use a relaxed SMARTS for a common pyranose sugar ring.
    glycone_smarts = "C1OC(C(O)C(O)C1O)"
    glycone_query = Chem.MolFromSmarts(glycone_smarts)
    if glycone_query is None:
        return None, None
    glycone_matches = mol.GetSubstructMatches(glycone_query, useChirality=False)
    if not glycone_matches:
        return False, "No glycone (sugar) ring found in the molecule."
    
    # --- Step 3: Look for a candidate central carbon that satisfies connectivity --
    # For each sulfonated oxime match we consider its central carbon (first atom in our oxime SMARTS).
    reasons = []
    for match in oxime_matches:
        # match[0] is the central carbon (the one with the C=N bond)
        central_c_idx = match[0]
        central_c = mol.GetAtomWithIdx(central_c_idx)
        
        # Check that the central carbon has at least three substituents (allows a side-chain)
        if central_c.GetDegree() < 3:
            reasons.append("Central carbon in sulfonated oxime has fewer than 3 substituents.")
            continue
        
        # Look for a neighbor sulfur (atomic number 16) that is not part of this oxime group.
        # In our oxime SMARTS the sulfonyl S is match[3] (if the match length is 4).
        candidate_sulfur = None
        for nbr in central_c.GetNeighbors():
            if nbr.GetAtomicNum() == 16:
                # if this sulfur is the one already in the oxime match, skip it
                if len(match) >= 4 and nbr.GetIdx() == match[3]:
                    continue
                candidate_sulfur = nbr
                break
        
        if candidate_sulfur is None:
            reasons.append("Central carbon is not attached to a sulfur atom (glycoside S) outside the oxime group.")
            continue
        
        # Now check that this candidate sulfur is connected (by at least one neighboring atom)
        # to a sugar (glycone) ring. We iterate over neighbors of the sulfur (skipping the central carbon).
        glycone_found = False
        for s_nbr in candidate_sulfur.GetNeighbors():
            if s_nbr.GetIdx() == central_c_idx:
                continue
            # Check if this neighbor atom belongs to any glycone match.
            for gly_match in glycone_matches:
                if s_nbr.GetIdx() in gly_match:
                    glycone_found = True
                    break
            if glycone_found:
                break
        
        if not glycone_found:
            reasons.append("Sulfur attached to central carbon is not connected to a glycone (sugar) ring.")
            continue
        
        # If we reach here for one of the matches, we consider the criteria satisfied.
        return True, ("Molecule contains a sulfonated oxime group on a central carbon (with a side chain) "
                      "and a sulfur connecting a pyranose (glycone) ring.")
    
    # If none of the oxime matches produce a valid candidate, return failure.
    full_reason = " ; ".join(reasons) if reasons else "No valid glucosinolate substructure found."
    return False, full_reason


# Example usage and testing:
if __name__ == "__main__":
    # Test with one of the provided glucosinolate SMILES:
    tests = [
        ("S(C1OC(C(O)C(O)C1O)CO\\C(=N\\OS(O)(=O)=O)\\CC/C=C/S(=O)C",
         "4-(Methylsulfinyl)but-3-enylglucosinolate"),
        ("S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CCC2=CC(=CC=C2)OC",
         "Phenethyl glucosinolate"),
        ("[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS([O-])(=O)=O)/CC/C=C/S(=O)C",
         "Another example")
    ]
    for smi, name in tests:
        result, reason = is_glucosinolate(smi)
        print(f"Testing {name}:\n SMILES: {smi}\n Result: {result}\n Reason: {reason}\n")