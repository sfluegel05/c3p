"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: 3-substituted propionyl-CoA(4-)

Definition:
 An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
 of any 3-substituted propionyl-CoA; major species at pH 7.3.

This implementation now checks:
 1. That the net formal charge is exactly –4.
 2. That the adenine moiety (the CoA headgroup) is present.
 3. That the characteristic CoA connector (pantetheine fragment) is found.
      A modified SMARTS is used: "SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP"
      (the use of "C(O)" rather than a chiral substructure makes the match more tolerant)
 4. That a thioester linkage is present – using an explicit "C(=O)S" pattern.
 5. That this thioester linkage’s sulfur is directly part of (or adjacent to) the CoA connector.
 6. That the acyl chain attached to the thioester carbonyl (and not part of the CoA scaffold)
    is composed solely of carbon atoms and is linear.
If any test fails, the molecule is not classified as a member of the 3‐substituted propionyl-CoA(4-) class.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule matches the class of 3-substituted propionyl-CoA(4-) (an acyl-CoA(4-) oxoanion)
    based on its SMILES string.
    
    Expected criteria:
      - Net formal charge exactly -4.
      - The characteristic adenine moiety (nucleobase) is present.
      - A composite CoA connector (pantetheine unit) is present.
          SMARTS used: "SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP"
      - A thioester linkage is present, defined as "C(=O)S", whose sulfur is part of (or directly adjacent
        to) the CoA connector.
      - The acyl chain attached at the thioester carbonyl (and not part of the CoA scaffold) is composed solely 
        of carbon atoms and is linear (no branching).
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        (bool, str): Tuple (True, reason) if the molecule is classified as 3‐substituted propionyl-CoA(4-);
                     otherwise (False, reason) explaining why.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Check net formal charge ---
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return False, f"Net formal charge is {net_charge} (expected -4)"
    
    # --- 2. Check for the adenine moiety using a fused purine SMARTS ---
    adenine_smarts = "n1cnc2cncnc12"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine):
        return False, "Adenine moiety not found"
    
    # --- 3. Check for the CoA connector (pantetheine) fragment ---
    # We use a slightly relaxed SMARTS so that chiral annotation does not exclude the match.
    coa_connector_smarts = "SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP"
    coa_connector = Chem.MolFromSmarts(coa_connector_smarts)
    coa_matches = mol.GetSubstructMatches(coa_connector)
    if not coa_matches:
        return False, "CoA connector (pantetheine fragment) not found"
    
    # Collect all atom indices from all CoA connector matches for cross-checking.
    # In addition, we will also allow atoms which are directly bonded to a connector atom.
    connector_atom_set = set()
    for match in coa_matches:
        connector_atom_set.update(match)
    # Also populate neighbors of the connector atoms (to allow slight mismatches)
    connector_neighbors = set()
    for idx in connector_atom_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            connector_neighbors.add(nbr.GetIdx())
    # We combine these into a set of acceptable indices for the thioester sulfur.
    acceptable_connector_indices = connector_atom_set.union(connector_neighbors)

    # --- 4. Check for the thioester linkage ---
    # Use a SMARTS that explicitly captures a carbonyl S connectivity
    thioester_smarts = "C(=O)S"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    ts_matches = mol.GetSubstructMatches(thioester)
    if not ts_matches:
        return False, "Thioester linkage (C(=O)S) not found"
    
    # Ensure that at least one thioester match has its sulfur atom connected to the connector.
    valid_ts = []
    for match in ts_matches:
        # In the match, index 0 is the carbonyl carbon, index 1 is the sulfur.
        carbonyl_idx, sulfur_idx = match[0], match[1]
        if sulfur_idx in acceptable_connector_indices:
            valid_ts.append(match)
    if not valid_ts:
        return False, "No thioester linkage found whose sulfur is part of or adjacent to the CoA connector"
    
    # --- 5. Validate the acyl chain attached at the thioester carbonyl ---
    # The acyl chain is defined as the carbon chain attached to the carbonyl carbon (other than the S).
    def extract_linear_chain_length(carbonyl_atom, exclude_idx):
        # Identify neighbor(s) of carbonyl that are carbons but not the excluded S atom.
        neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() 
                     if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != exclude_idx]
        if not neighbors:
            return None  # No acyl branch found.
        # We assume a single acyl chain connection.
        start = neighbors[0]
        length = 1  # count the first carbon in the acyl chain
        prev = carbonyl_atom
        current = start
        while True:
            # Look for a unique carbon neighbor not coming from the previous atom.
            next_neighbors = [nbr for nbr in current.GetNeighbors() 
                              if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev.GetIdx()]
            if len(next_neighbors) == 0:
                break  # reached the chain end
            elif len(next_neighbors) == 1:
                length += 1
                prev, current = current, next_neighbors[0]
            else:
                # Branching detected; invalid acyl chain.
                return None
        return length

    acyl_chain_valid = False
    acyl_chain_length = None
    for ts in valid_ts:
        carbonyl_idx, sulfur_idx = ts[0], ts[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Get the neighbor carbons of the carbonyl (excluding the thioester S)
        acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() 
                          if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != sulfur_idx]
        if not acyl_neighbors:
            continue  # no acyl chain attached
        chain_len = extract_linear_chain_length(carbonyl_atom, sulfur_idx)
        if chain_len is None:
            continue  # invalid (branched or non-carbon) chain
        acyl_chain_valid = True
        acyl_chain_length = chain_len
        break
    if not acyl_chain_valid:
        return False, "Acyl chain is branched or contains non-carbon atoms"
    
    # Optionally enforce a minimal chain length (here, we expect at least 3 carbons)
    if acyl_chain_length < 3:
        return False, f"Acyl chain is too short (length {acyl_chain_length}); expected a fatty acyl chain"

    return True, ("Molecule matches structural criteria for 3-substituted propionyl-CoA(4-); "
                  f"acyl chain length (excluding carbonyl) = {acyl_chain_length}")

# Example usage:
if __name__ == "__main__":
    # Example SMILES for tetracosanoyl-CoA(4-)
    smiles_example = ("CCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)"
                      "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)"
                      "OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)"
                      "n1cnc2c(N)ncnc12")
    result, reason = is_3_substituted_propionyl_CoA_4__(smiles_example)
    print(result, reason)