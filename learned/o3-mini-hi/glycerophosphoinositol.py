"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: Glycerophosphoinositol

Definition:
  Any glycerophospholipid having the polar alcohol inositol esterified to the phosphate group 
  at the sn-3 position of the glycerol backbone.

The algorithm (improved):
  1. Parse the SMILES string.
  2. Search for a glycerol backbone using an approximate fragment ("OCC(O)CO").
  3. Search for an inositol head‐group using an approximate fully hydroxylated cyclohexane pattern.
  4. Loop over phosphorus atoms in the molecule:
      • Reject any phosphorus that is “polysubstituted” (i.e. any oxygen attached to P that also 
        bonds to another phosphorus, to avoid polyphosphate species).
      • For each remaining phosphorus, check whether at least one oxygen neighbor is next to a 
        glycerol atom and at least one is next to an inositol atom.
  5. Return True only if a “bridging” phosphate is found.
  
Note: The SMARTS/SMILES queries for glycerol and inositol are only approximate. Additional 
      constraints (such as ensuring the phosphate is attached at the sn-3 position) would require
      a more sophisticated pattern.
      
If any step fails (e.g. invalid SMILES or required substructure not found), the function returns 
(False, reason). If it is unable to confidently classify then (None, None) may be returned.
"""

from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is glycerophosphoinositol based on its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as glycerophosphoinositol, False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Search for a glycerol backbone.
    # We use an approximate fragment for glycerol: HO-CH2-CHOH-CH2-OH, as "OCC(O)CO".
    glycerol_query = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_query is None:
        return False, "Error creating glycerol query"
    glycerol_matches = mol.GetSubstructMatches(glycerol_query)
    if not glycerol_matches:
        return False, "No glycerol backbone fragment found"
    # Record all atom indices that participate in any glycerol match.
    glycerol_atoms = {idx for match in glycerol_matches for idx in match}
    
    # 3. Search for an inositol head‐group.
    # Use an approximate pattern: fully hydroxylated cyclohexane ("C1(C(C(C(C(C1)O)O)O)O)O").
    inositol_query = Chem.MolFromSmiles("C1(C(C(C(C(C1)O)O)O)O)O")
    if inositol_query is None:
        return False, "Error creating inositol query"
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "No inositol head-group found"
    inositol_atoms = {idx for match in inositol_matches for idx in match}
    
    # 4. Now search for a bridging phosphate group.
    # Loop through any phosphorus (atomic number 15) in the molecule.
    phosphate_bridge_found = False
    phosphorus_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        phosphorus_found = True
        # Gather oxygen neighbors.
        oxy_neighbors = []
        skip_due_to_polyphosphate = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen
                oxy_neighbors.append(nbr)
                # If an oxygen is connected to any other phosphorus, skip this candidate.
                for nn in nbr.GetNeighbors():
                    if nn.GetAtomicNum() == 15 and nn.GetIdx() != atom.GetIdx():
                        skip_due_to_polyphosphate = True
                        break
            if skip_due_to_polyphosphate:
                break
        if skip_due_to_polyphosphate:
            continue  # This phosphorus is part of a polyphosphate (or extra phosphate substitution).
        
        # Require at least two oxygen neighbors.
        if len(oxy_neighbors) < 2:
            continue
        
        # Check if the phosphate connects to both glycerol and inositol fragments.
        connects_to_glycerol = False
        connects_to_inositol = False
        for oxy in oxy_neighbors:
            for nbr in oxy.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetIdx() in glycerol_atoms:
                    connects_to_glycerol = True
                if nbr.GetIdx() in inositol_atoms:
                    connects_to_inositol = True
            # Early exit if both connections are seen.
            if connects_to_glycerol and connects_to_inositol:
                break
        
        if connects_to_glycerol and connects_to_inositol:
            phosphate_bridge_found = True
            break
    
    if not phosphorus_found:
        return False, "No phosphorus atom found in molecule"
    if not phosphate_bridge_found:
        return False, "No phosphate group found that cleanly bridges glycerol and inositol fragments"
    
    return True, "Molecule contains a glycerol backbone with a phosphate group bridging to an inositol head-group"

# Example usage for testing:
if __name__ == "__main__":
    # Test example (PHOOA-PI)
    test_smiles = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCC(O)/C=C/C=O)COC(=O)CCCCCCCCCCCCCCC)(O)=O"
    result, explanation = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Explanation:", explanation)