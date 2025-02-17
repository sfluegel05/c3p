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

    Strategy:
      1. Look for a sulfonated oxime group. Our modified SMARTS pattern now expects:
         a carbon double bonded to an N that is connected to an O which in turn
         is bonded to an S and two oxygens in a sulfonyl moiety.
         SMARTS: "C(=N[O]S(=O)(=O)[O-]?)" (the charge on the terminal O is optional).
      2. Verify that the central carbon (the one in the oxime) has at least three substituents.
      3. From the central carbon, find a neighbor sulfur atom that is not the sulfonyl S
         (i.e. not the one coming from the oxime match).
      4. For that candidate sulfur, search its neighbors (ignoring the central carbon)
         to see if one of them is part of a sugar (glycone) ring.
         We use a relaxed SMARTS for a pyranose ring: "C1OC(C(O)C(O)C1O)"
         and ignore chirality.
      5. If a match is found, return True with an explanation.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as glucosinolate, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- Step 1: Find the sulfonated oxime group ---
    # Our modified SMARTS pattern for the sulfonated oxime fragment.
    oxime_smarts = "C(=N[O]S(=O)(=O)[O-]?)"
    oxime_query = Chem.MolFromSmarts(oxime_smarts)
    if oxime_query is None:
        return None, None  # unexpected error with SMARTS pattern
    
    # We ignore chirality in the substructure search.
    oxime_matches = mol.GetSubstructMatches(oxime_query, useChirality=False)
    if not oxime_matches:
        return False, "No sulfonated oxime group found."
    
    # --- Step 2: Find the glycone (sugar) ring ---
    # Use a relaxed SMARTS pattern for a common pyranose sugar ring.
    glycone_smarts = "C1OC(C(O)C(O)C1O)"
    glycone_query = Chem.MolFromSmarts(glycone_smarts)
    if glycone_query is None:
        return None, None
    
    glycone_matches = mol.GetSubstructMatches(glycone_query, useChirality=False)
    if not glycone_matches:
        return False, "No glycone (sugar) ring found in the molecule."
    
    # --- Step 3: For each sulfonated oxime match, check connectivity ---
    # In our oxime SMARTS, we expect the following indexing:
    # match[0] = central carbon (C), match[1] = nitrogen (N), match[2] = oxygen (O) of the oxime,
    # match[3] = sulfonyl sulfur (S) [if matched]. This assumes the pattern is matched as written.
    reasons = []
    for match in oxime_matches:
        # Ensure the match has at least 3 atoms (our pattern should yield 4 atoms if the sulfonyl S is included)
        if len(match) < 3:
            reasons.append("Oxime match did not yield expected number of atoms.")
            continue
        
        central_c_idx = match[0]
        central_c = mol.GetAtomWithIdx(central_c_idx)
        
        # The central C should have degree at least 3 to have a side-chain.
        if central_c.GetDegree() < 3:
            reasons.append("Central carbon in sulfonated oxime has fewer than 3 substituents.")
            continue
        
        # Find a neighbor sulfur atom attached to the central carbon that is NOT the sulfonyl S from the oxime.
        candidate_sulfur = None
        for nbr in central_c.GetNeighbors():
            if nbr.GetAtomicNum() == 16:
                # If our oxime match included a sulfonyl S (match[3] if available) then skip that atom.
                if len(match) >= 4 and nbr.GetIdx() == match[3]:
                    continue
                candidate_sulfur = nbr
                break
        
        if candidate_sulfur is None:
            reasons.append("Central carbon is not attached to a sulfur atom (glycosidic S) outside the oxime group.")
            continue
        
        # Check that candidate sulfur is connected to a glycone (sugar) ring.
        glycone_found = False
        for s_nbr in candidate_sulfur.GetNeighbors():
            # Skip re-checking connection to the central carbon.
            if s_nbr.GetIdx() == central_c_idx:
                continue
            # If this neighbor is part of any glycone ring match, count it.
            for gly_match in glycone_matches:
                if s_nbr.GetIdx() in gly_match:
                    glycone_found = True
                    break
            if glycone_found:
                break
        
        if not glycone_found:
            reasons.append("The candidate sulfur attached to the central carbon is not connected to a glycone (sugar) ring.")
            continue
        
        # If all criteria are met for one match, we have a valid glucosinolate.
        return True, ("Molecule contains a sulfonated oxime fragment on a central carbon with sufficient substituents, "
                      "and this central carbon is bonded to a sulfur which connects to a glycone (sugar) ring.")
    
    # If none of the matches produce a valid candidate, return False with the aggregated reasons.
    full_reason = " ; ".join(reasons) if reasons else "No valid glucosinolate substructure found."
    return False, full_reason


# Example usage and testing:
if __name__ == "__main__":
    # List of example SMILES for glucosinolates.
    test_smiles = [
        "S(C1OC(C(O)C(O)C1O)COC(=O)/C=C/C2=CC(OC)=C(O)C(OC)=C2)\\C(=N\\OS(O)(=O)=O)\\CC/C=C/S(=O)C",  # example 1
        "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CC/C=C/S(=O)C",  # example 2
        "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS([O-])(=O)=O)/CC2=CC(=CC=C2)OC",  # example 3
        "OC[C@H]1O[C@@H](S\\C(CCCCOC(=O)C2=CC=CC=C2)=N\\OS([O-])(=O)=O)[C@H](O)[C@@H](O)[C@@H]1O",  # example 4
        "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CCC2=CC(=CC=C2)",  # Phenethyl glucosinolate
    ]
    
    for i, smi in enumerate(test_smiles, start=1):
        result, reason = is_glucosinolate(smi)
        print(f"Test {i}:\n SMILES: {smi}\n Result: {result}\n Reason: {reason}\n")