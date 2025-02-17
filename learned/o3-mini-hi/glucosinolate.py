"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: Glucosinolate
Definition: Water‐soluble anionic substituted thioglucosides.
Glucosinolates have a central C atom bonded via an S atom to a glycone (sugar) group and via an N atom 
to a sulfonated oxime group; this central carbon also carries a side‐group.
(The anti configuration is not explicitly verified.)
"""

from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    Strategy:
      1. Identify a sulfonated oxime fragment using a SMARTS pattern. This pattern
         looks for a central carbon double-bonded to a nitrogen that is in turn bonded
         to an oxygen and sulfur (S(=O)(=O)[O–] optionally) i.e. "C(=N[O]S(=O)(=O)[O-]?)".
      2. Identify a glycone (sugar) fragment using a SMARTS pattern for a pyranose-like ring.
         (Here we use a common motif "C1OC(C(O)C(O)C1O)".)
      3. For each oxime match, check that its central carbon has at least 3 substituents.
      4. Among those substituents, find a sulfur (other than the sulfur that is part of the oxime fragment).
      5. Check if that candidate sulfur is connected to any atom that is part of a glycone match.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a glucosinolate, False otherwise.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- Step 1: Search for the sulfonated oxime fragment ---
    # This SMARTS looks for a fragment with a carbon double-bonded to an N that is bonded to an oxygen,
    # which is in turn attached to a sulfonyl group S(=O)(=O)[O-] (the negative charge is optional).
    oxime_smarts = "C(=N[O]S(=O)(=O)[O-]?)"
    oxime_query = Chem.MolFromSmarts(oxime_smarts)
    if oxime_query is None:
        return None, None  # error parsing SMARTS
    
    # Use useChirality=False to be more general.
    oxime_matches = mol.GetSubstructMatches(oxime_query, useChirality=False)
    if not oxime_matches:
        return False, "No sulfonated oxime fragment found."
    
    # --- Step 2: Search for the glycone (sugar) ring ---
    # We use a SMARTS pattern that is common for a pyranose ring.
    # Note: many glucosinolate sugar parts are derived from glucose.
    glycone_smarts = "C1OC(C(O)C(O)C1O)"
    glycone_query = Chem.MolFromSmarts(glycone_smarts)
    if glycone_query is None:
        return None, None  # error parsing SMARTS
    glycone_matches = mol.GetSubstructMatches(glycone_query, useChirality=False)
    if not glycone_matches:
        return False, "No glycone (sugar) ring found in the molecule."
    # Gather all atom indices that are part of any glycone match.
    glycone_atoms = set()
    for match in glycone_matches:
        glycone_atoms.update(match)
    
    # --- Step 3: For each oxime match, check for the connectivity criteria ---
    reasons = []  # collect failure reasons from various matches
    for match in oxime_matches:
        # Our oxime SMARTS is expected to match 4 atoms:
        # match[0]: central carbon (C)
        # match[1]: imine nitrogen (N)
        # match[2]: oxime oxygen (O)
        # match[3]: sulfonyl sulfur (S) within the oxime fragment.
        if len(match) < 4:
            reasons.append("Oxime match did not return the expected four atoms.")
            continue
        central_c_idx = match[0]
        central_c = mol.GetAtomWithIdx(central_c_idx)
        # The central carbon should have at least 3 substituents (one for the oxime nitrogen, one for the candidate S, and one side-chain)
        if central_c.GetDegree() < 3:
            reasons.append("Central carbon in the oxime fragment has fewer than 3 substituents.")
            continue
        
        # --- Step 4: Look for a candidate sulfur attached to the central carbon
        candidate_sulfur = None
        for nbr in central_c.GetNeighbors():
            if nbr.GetAtomicNum() == 16:  # atomic number for sulfur is 16
                # Exclude the sulfur that is part of the oxime fragment (match[3])
                if nbr.GetIdx() == match[3]:
                    continue
                candidate_sulfur = nbr
                break
        if candidate_sulfur is None:
            reasons.append("Central carbon is not attached to a second sulfur atom (needed for glycone linkage).")
            continue
        
        # --- Step 5: Check if the candidate sulfur is connected to a glycone (sugar) ring.
        glycone_found = False
        # Look at all neighbors of the candidate sulfur (except the central carbon):
        for nbr in candidate_sulfur.GetNeighbors():
            if nbr.GetIdx() == central_c_idx:
                continue
            if nbr.GetIdx() in glycone_atoms:
                glycone_found = True
                break
                
        if not glycone_found:
            reasons.append("The candidate sulfur attached to the central carbon is not connected to a glycone (sugar) ring.")
            continue
        
        # If we reach here, the molecule meets the key criteria for a glucosinolate.
        return True, ("Molecule contains a sulfonated oxime fragment (with pattern C(=N[O]S(=O)(=O)[O-]?) and a central carbon "
                      "that is attached to a candidate sulfur, which in turn is connected to a glycone (sugar) ring.")
    
    # If none of the alternative oxime matches lead to a proper candidate, return the reasons collected.
    full_reason = " ; ".join(reasons) if reasons else "No valid glucosinolate substructure found."
    return False, full_reason


# (Optional) Example usage:
if __name__ == "__main__":
    # A few test SMILES from the provided list.
    test_smiles = [
        "S(C1OC(C(O)C(O)C1O)COC(=O)/C=C/C2=CC(OC)=C(O)C(OC)=C2)\\C(=N\\OS(O)(=O)=O)\\CC/C=C/S(=O)C",
        "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CC/C=C/S(=O)C",
        "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS([O-])(=O)=O)/CC2=CC(=CC=C2)OC",
        "OC[C@H]1O[C@@H](S\\C(CCCCOC(=O)C2=CC=CC=C2)=N\\OS([O-])(=O)=O)[C@H](O)[C@@H](O)[C@@H]1O",
        "S(C1OC(C(O)C(O)C1O)CO)\\C(=N\\OS(O)(=O)=O)\\CCC2=CC(=CC=C2)"
    ]
    
    for i, smi in enumerate(test_smiles, start=1):
        result, reason = is_glucosinolate(smi)
        print(f"Test {i}:\n SMILES: {smi}\n Result: {result}\n Reason: {reason}\n")