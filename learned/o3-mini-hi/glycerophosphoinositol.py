"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: Glycerophosphoinositol

Definition:
  Any glycerophospholipid having the polar alcohol inositol esterified to the phosphate group 
  at the sn-3 position of the glycerol backbone.

The algorithm:
  1. Parse the SMILES & check that we have a valid molecule.
  2. Search for a glycerol backbone by matching a simple pattern (approximated as "OCC(O)CO").
  3. Search for an inositol head-group using an approximated fully hydroxylated cyclohexane ring.
  4. Loop over phosphorus atoms (P) in the molecule:
       For each phosphorus, get its oxygen neighbors.
       For each oxygen neighbor, check its other neighbor(s) (ignoring phosphorus).
       If at least one such oxygen is attached to any atom that is part of the glycerol group 
         and if another oxygen (or a different neighbor) is attached to any atom that is part of 
         the inositol head-group, then we say that the phosphate bridges glycerol and inositol.
       
  If all criteria are met, the molecule is classified as glycerophosphoinositol.
"""

from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is glycerophosphoinositol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as glycerophosphoinositol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for glycerol backbone.
    # This is approximated by a fragment: HO-CH2-CHOH-CH2-OH, represented as "OCC(O)CO".
    glycerol_query = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_query is None:
        return False, "Error creating glycerol query"
    glycerol_matches = mol.GetSubstructMatches(glycerol_query)
    if not glycerol_matches:
        return False, "No glycerol backbone fragment found"
    # Record all atom indices involved in any glycerol match.
    glycerol_atoms = {idx for match in glycerol_matches for idx in match}
    
    # 2. Look for an inositol head-group.
    # Using an approximate pattern: fully hydroxylated cyclohexane.
    # Chirality is not enforced for the query so that different representations match.
    inositol_query = Chem.MolFromSmiles("C1(C(C(C(C(C1)O)O)O)O)O")
    if inositol_query is None:
        return False, "Error creating inositol query"
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "No inositol head-group found"
    inositol_atoms = {idx for match in inositol_matches for idx in match}
    
    # 3. Now search for a bridging phosphate group.
    # We loop through every phosphorus atom (atomic number 15)
    phosphate_bridge_found = False
    phosphorus_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        phosphorus_found = True

        # Check all oxygen neighbors of the phosphorus.
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen
                oxy_neighbors.append(nbr)
                
        # We need at least two oxygens attached to phosphorus
        if len(oxy_neighbors) < 2:
            continue
        
        # For each oxygen neighbor, check if it connects (via its other neighbor) to glycerol or inositol.
        connects_to_glycerol = False
        connects_to_inositol = False
        for oxy in oxy_neighbors:
            # For each oxygen, get its neighbors (excluding the phosphorus).
            for nbr in oxy.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                # If the neighbor atom is in the glycerol fragment, note that.
                if nbr.GetIdx() in glycerol_atoms:
                    connects_to_glycerol = True
                # If the neighbor atom is in the inositol fragment, note that.
                if nbr.GetIdx() in inositol_atoms:
                    connects_to_inositol = True
            # Early exit if both connections are seen.
            if connects_to_glycerol and connects_to_inositol:
                break
        
        if connects_to_glycerol and connects_to_inositol:
            phosphate_bridge_found = True
            break

    if not phosphorus_found:
        return False, "No phosphorus (phosphate) atom found in molecule"
    if not phosphate_bridge_found:
        return False, "No phosphate group found that bridges glycerol and inositol fragments"
    
    return True, "Molecule contains a glycerol backbone with a phosphate group bridging to an inositol head-group"

# Example usage for testing:
if __name__ == "__main__":
    # Example test (PHOOA-PI)
    test_smiles = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCC(O)/C=C/C=O)COC(=O)CCCCCCCCCCCCCCC)(O)=O"
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)