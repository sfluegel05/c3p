"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: Glycerophosphoinositol
Definition: Any glycerophospholipid having the polar alcohol inositol esterified 
to the phosphate group at the sn-3 position of the glycerol backbone.

This implementation uses improved substructure queries and stricter requirements:
  1. The molecule must contain a glycerol fragment (approximated by a glycerol SMARTS).
  2. The molecule must contain an inositol head-group (approximated by a fully hydroxylated cyclohexane).
  3. At least one phosphorus atom (P) must have a double-bonded oxygen (P=O) and exactly two single-bonded oxygen neighbors.
     Among these oxygens, one must be “clean” (attached only to that phosphorus) and connected to a carbon that lies in the glycerol fragment,
     and the other must be “clean” and connected to a carbon that lies in the inositol fragment.
     
If any of these conditions fail, the molecule is not classified as glycerophosphoinositol.
"""

from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    
    The algorithm checks that the molecule:
      1. Contains a glycerol backbone fragment (approximated with SMARTS "OCC(O)CO").
      2. Contains an inositol head-group (approximated by "C1(C(C(C(C(C1O)O)O)O)O)O").
      3. Contains at least one phosphate (P) atom where:
            - There is one oxygen neighbor with a double bond to P (the P=O) (or with bond order >1).
            - There are exactly two additional “ester-type” oxygen neighbors (with a single bond to P)
              that are “clean” (i.e. each is attached only to that phosphorus) 
              and where one oxygen connects (through its other neighbor) to the glycerol fragment 
              and the other connects to the inositol fragment.
    
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
    
    # Query for glycerol: approximated as HOCH2-CHOH-CH2OH. (SMARTS: OCC(O)CO)
    glycerol_query = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_query is None:
        return False, "Error creating glycerol query"
    if not mol.HasSubstructMatch(glycerol_query):
        return False, "No glycerol backbone fragment found"
    
    # Query for inositol: approximated as fully hydroxylated cyclohexane.
    inositol_query = Chem.MolFromSmiles("C1(C(C(C(C(C1O)O)O)O)O)O")
    if inositol_query is None:
        return False, "Error creating inositol query"
    if not mol.HasSubstructMatch(inositol_query):
        return False, "No inositol head-group found"
    
    # Collect atom indices of the glycerol and inositol matches.
    glycerol_atoms = set()
    for match in mol.GetSubstructMatches(glycerol_query):
        glycerol_atoms.update(match)
    inositol_atoms = set()
    for match in mol.GetSubstructMatches(inositol_query):
        inositol_atoms.update(match)
    
    # Now, loop over phosphorus atoms.
    phosphate_bridge_found = False
    phosphorus_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        phosphorus_found = True
        
        # Collect oxygen neighbors and record bond orders.
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # Save tuple of (oxygen atom, bond order)
                order = bond.GetBondTypeAsDouble()
                oxy_neighbors.append((nbr, order))
        
        # We expect one double bond (P=O) and exactly two single bonds to oxygen that serve as ester linkages.
        # (Other oxygens, e.g. -OH, might be present but we focus on the two linking oxygens.)
        double_bond_count = sum(1 for oxy, oorder in oxy_neighbors if oorder > 1.5)
        single_ester_candidates = [oxy for oxy, oorder in oxy_neighbors if abs(oorder - 1.0) < 0.1]
        
        # Require at least one double-bonded oxygen and at least two single-bonded oxygen neighbors.
        if double_bond_count < 1 or len(single_ester_candidates) < 2:
            continue
        
        # For each single-bond oxygen candidate, check if it is "clean": i.e. it is attached to only one phosphorus.
        clean_ester_oxygens = []
        for oxy in single_ester_candidates:
            p_count = 0
            for nbr in oxy.GetNeighbors():
                if nbr.GetAtomicNum() == 15:
                    p_count += 1
            if p_count == 1:
                clean_ester_oxygens.append(oxy)
        
        if len(clean_ester_oxygens) < 2:
            continue  # Possibly a polyphosphate situation; skip this P.
        
        # Now, determine if among these two clean ester oxygens one connects (through its other neighbor)
        # to an atom in the glycerol fragment and the other to an atom in the inositol fragment.
        connects_glycerol = False
        connects_inositol = False
        for oxy in clean_ester_oxygens:
            for nbr in oxy.GetNeighbors():
                # Skip the phosphorus itself.
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                # We want the oxygen to be attached to a carbon (or heteroatom) that belongs to one of our fragments.
                if nbr.GetAtomicNum() == 6:  # carbon atoms in organic fragments
                    if nbr.GetIdx() in glycerol_atoms:
                        connects_glycerol = True
                    if nbr.GetIdx() in inositol_atoms:
                        connects_inositol = True
                else:
                    # Sometimes the connection might be via another oxygen,
                    # check one level deeper.
                    for nn in nbr.GetNeighbors():
                        if nn.GetIdx() == oxy.GetIdx() or nn.GetIdx() == atom.GetIdx():
                            continue
                        if nn.GetIdx() in glycerol_atoms:
                            connects_glycerol = True
                        if nn.GetIdx() in inositol_atoms:
                            connects_inositol = True
            if connects_glycerol and connects_inositol:
                break
        
        if connects_glycerol and connects_inositol:
            phosphate_bridge_found = True
            break
    
    if not phosphorus_found:
        return False, "No phosphorus (phosphate) atom found in molecule"
    if not phosphate_bridge_found:
        return False, "No phosphate group found that cleanly bridges the glycerol backbone and inositol head-group"
    
    return True, "Molecule contains a glycerol backbone with a phosphate group that directly connects to an inositol head-group"

# Example usage for testing:
if __name__ == "__main__":
    # Test with one example of a phosphatidylinositol (PI) structure.
    test_smiles = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCC(O)/C=C/C=O)COC(=O)CCCCCCCCCCCCCCC)(O)=O"
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)