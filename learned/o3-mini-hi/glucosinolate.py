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
      1. Look for a sulfonated oxime fragment. In our SMARTS, we expect a pattern:
         central carbon double-bonded to nitrogen, which is bonded to an oxygen that is linked
         to a sulfonyl group; i.e. "C(=N[O]S(=O)(=O)[O-]?)". This pattern is designed to capture the
         oxime fragment while allowing an optional negative charge.
      2. For each oxime match we examine the identified central carbon (match[0]). This carbon should
         have at least three substituents (to allow a side-chain).
      3. Among its neighbors, find a sulfur atom that is not the sulfonyl sulfur already in the oxime match.
      4. For that candidate sulfur, check whether one of its neighbors (other than the central C)
         is part of a glycone (sugar) ring. We use a standard pyranose pattern: "C1OC(C(O)C(O)C1O)"
         (ignoring chirality).
      5. If all criteria are met, the molecule is classified as a glucosinolate.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a glucosinolate, False otherwise.
        str: Explanation for the classification decision.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # --- Step 1: Find the sulfonated oxime fragment ---
    # The pattern below is designed to detect a C=NO-S(=O)(=O)[O-]? fragment.
    oxime_smarts = "C(=N[O]S(=O)(=O)[O-]?)"
    oxime_query = Chem.MolFromSmarts(oxime_smarts)
    if oxime_query is None:
        return None, None  # error with SMARTS
    
    # ignore chirality since many examples use stereochemistry
    oxime_matches = mol.GetSubstructMatches(oxime_query, useChirality=False)
    if not oxime_matches:
        return False, "No sulfonated oxime fragment found."

    # --- Step 2: Find the glycone (sugar) ring ---
    # The following represents a common pyranose ring substructure.
    glycone_smarts = "C1OC(C(O)C(O)C1O)"
    glycone_query = Chem.MolFromSmarts(glycone_smarts)
    if glycone_query is None:
        return None, None
    glycone_matches = mol.GetSubstructMatches(glycone_query, useChirality=False)
    if not glycone_matches:
        return False, "No glycone (sugar) ring found in the molecule."

    # --- Step 3: For each oxime match check connectivity ---
    # We expect the oxime SMARTS to return atoms in this order:
    # match[0]: central carbon (C) that is double bonded to nitrogen.
    # match[1]: imine nitrogen (N)
    # match[2]: oxime oxygen (O)
    # match[3]: sulfonyl sulfur (S) (part of the oxime fragment)
    reasons = []
    for match in oxime_matches:
        # We expect at least 4 atoms from our oxime to know the sulfonyl S index.
        if len(match) < 4:
            reasons.append("Oxime match did not yield expected number of atoms.")
            continue

        central_c_idx = match[0]
        central_c = mol.GetAtomWithIdx(central_c_idx)
        # Require that the central carbon has at least three substituents (for side-group attachment)
        if central_c.GetDegree() < 3:
            reasons.append("Central carbon in the oxime fragment has fewer than 3 substituents.")
            continue

        # --- Step 4: Look for a candidate sulfur attached to the central C that is not the sulfonyl S.
        candidate_sulfur = None
        for nbr in central_c.GetNeighbors():
            if nbr.GetAtomicNum() == 16:  # sulfur atomic number
                # Exclude the sulfur that is part of the oxime match (our match[3])
                if nbr.GetIdx() == match[3]:
                    continue
                candidate_sulfur = nbr
                break
        if candidate_sulfur is None:
            reasons.append("Central carbon is not attached to a second sulfur atom (expected for glycone linkage).")
            continue

        # --- Step 5: Check that candidate sulfur is connected (directly) to a glycone (sugar) ring.
        glycone_found = False
        for cand_nbr in candidate_sulfur.GetNeighbors():
            # Skip the central carbon
            if cand_nbr.GetIdx() == central_c_idx:
                continue
            # Check if this neighbor atom is part of any glycone ring match.
            for gly_match in glycone_matches:
                if cand_nbr.GetIdx() in gly_match:
                    glycone_found = True
                    break
            if glycone_found:
                break

        if not glycone_found:
            reasons.append("The sulfur substituent off the central carbon is not connected to a glycone (sugar) ring.")
            continue

        # If one instance meets all criteria, we classify as glucosinolate.
        return True, ("Molecule contains a sulfonated oxime fragment (central C=NO-S(=O)(=O)[O-]?) with a central carbon that "
                      "has a candidate glycosidic sulfur connected to a sugar ring.")
    
    # If none of the oxime matches lead to a valid candidate, then return failure with reasons.
    full_reason = " ; ".join(reasons) if reasons else "No valid glucosinolate substructure found."
    return False, full_reason


# Example usage and testing (can be removed or expanded as needed):
if __name__ == "__main__":
    # A small selection of example SMILES from the provided list.
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