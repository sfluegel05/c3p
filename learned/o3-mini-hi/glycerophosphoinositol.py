"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: Glycerophosphoinositol

Definition: Any glycerophospholipid having the polar alcohol inositol esterified 
to the phosphate group at the sn-3 position of the glycerol backbone.

Improvement in this version:
  1. The molecule must contain a glycerol backbone fragment (approximated here by SMARTS "OCC(O)CO").
  2. The molecule must contain an inositol head-group (approximated by using the SMILES for a fully hydroxylated cyclohexane ring).
  3. We then search for a candidate phosphate atom (atomic number 15) that meets the following:
       • It has exactly 3 oxygen neighbors.
       • Exactly one oxygen is double bonded (P=O; bond order > 1.5).
       • The two remaining oxygen neighbors are “ester-type”; each should be “clean” (i.e. have exactly 2 neighbors, the phosphorus and one other atom).
       • One of those oxygen’s other neighbor belongs to a glycerol fragment while the other oxygen’s neighbor belongs to an inositol fragment.
       
If any of these conditions fail, the molecule is not classified as glycerophosphoinositol.
"""

from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES string into an RDKit Mol.
      2. Checks for the presence of a glycerol backbone, approximated by the SMARTS "OCC(O)CO".
      3. Checks for the presence of an inositol head-group, approximated by the SMILES for fully hydroxylated cyclohexane.
      4. Loops over candidate phosphorus atoms (P) in the molecule. For each phosphorus, it
         - Collects oxygen neighbors and requires exactly three.
         - Requires exactly one double-bonded oxygen (P=O) (bond order > 1.5).
         - Among the two remaining single-bonded oxygen neighbors, each must be “clean” (have exactly 2 neighbors).
         - Checks that one of these oxygens is attached (via its other neighbor) directly to a glycerol fragment,
           and the other is attached to an inositol fragment.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as glycerophosphoinositol, False otherwise.
        str: Reason for the classification decision.
    """
    # Convert SMILES to molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Query for glycerol backbone fragment: approximated as HOCH2-CHOH-CH2OH ("OCC(O)CO")
    glycerol_query = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_query is None:
        return False, "Error creating glycerol query"
    glycerol_matches = mol.GetSubstructMatches(glycerol_query)
    if not glycerol_matches:
        return False, "No glycerol backbone fragment found"
    # Record all atom indices in glycerol matches.
    glycerol_atoms = {idx for match in glycerol_matches for idx in match}
    
    # Query for inositol head-group: approximated by fully hydroxylated cyclohexane.
    # (This is a simplification; many protonation states are possible.)
    inositol_query = Chem.MolFromSmiles("C1(C(C(C(C(C1O)O)O)O)O)O")
    if inositol_query is None:
        return False, "Error creating inositol query"
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "No inositol head-group found"
    inositol_atoms = {idx for match in inositol_matches for idx in match}
    
    # Now, search for a bridging phosphate.
    # We require a phosphorus atom (atomic number 15) that has exactly 3 oxygen neighbors:
    # one double-bonded (P=O) and two single-bonded "ester" oxygens.
    phosphate_bridge_found = False
    phosphorus_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        phosphorus_found = True
        
        # Get all oxygen neighbors of the phosphorus along with bond order.
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None:
                continue
            order = bond.GetBondTypeAsDouble()
            oxy_neighbors.append((nbr, order))
        
        # Require exactly three oxygen neighbors.
        if len(oxy_neighbors) != 3:
            continue
        
        # Count how many are double bonded.
        dbl_oxy = [oxy for oxy, oorder in oxy_neighbors if oorder > 1.5]
        if len(dbl_oxy) != 1:
            continue
        
        # The remaining two must be single-bonded ester oxygens.
        ester_candidates = [oxy for oxy, oorder in oxy_neighbors if abs(oorder - 1.0) < 0.1]
        if len(ester_candidates) != 2:
            continue
        
        # For each ester oxygen candidate, check it is “clean” (only bonded to phosphorus and one other atom).
        clean_ester_oxygens = []
        for oxy in ester_candidates:
            # Count neighbors not counting the phosphorus.
            nbrs = [n for n in oxy.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
            if len(nbrs) == 1:
                clean_ester_oxygens.append((oxy, nbrs[0]))
        if len(clean_ester_oxygens) != 2:
            continue
        
        # Now check connectivity:
        # We require that one of the ester oxygens connects directly to an atom in the glycerol fragment
        # and the other connects directly to an atom in the inositol fragment.
        connects_glycerol = False
        connects_inositol = False
        for oxy, ext_atom in clean_ester_oxygens:
            # We expect the connecting atom to be carbon.
            if ext_atom.GetAtomicNum() != 6 and ext_atom.GetAtomicNum() != 8:
                continue
            # Direct connection: if the external atom is part of the glycerol or inositol match.
            if ext_atom.GetIdx() in glycerol_atoms:
                connects_glycerol = True
            if ext_atom.GetIdx() in inositol_atoms:
                connects_inositol = True
        if connects_glycerol and connects_inositol:
            phosphate_bridge_found = True
            break

    if not phosphorus_found:
        return False, "No phosphorus (phosphate) atom found in molecule"
    if not phosphate_bridge_found:
        return False, "No phosphate group found that cleanly bridges glycerol and inositol fragments"
    
    return True, "Molecule contains a glycerol backbone with a phosphate group linking to an inositol head-group"

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCC(O)/C=C/C=O)COC(=O)CCCCCCCCCCCCCCC)(O)=O"
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)