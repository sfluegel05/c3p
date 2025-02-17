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
 3. That the characteristic CoA connector (pantetheine fragment) is found – here using a more tolerant SMARTS:
      "SCCNC(=O)CCNC(=O)[CX4](O)C(C)(C)COP"
 4. That a thioester linkage (pattern "[#6](=O)S") is present and that the sulfur in that linkage 
    is embedded within one of the CoA connector matches.
 5. That the acyl chain attached to the thioester carbonyl (that is not part of the CoA scaffold) 
    is made entirely of carbon atoms and is linear (no branching).
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
          SMARTS used: "SCCNC(=O)CCNC(=O)[CX4](O)C(C)(C)COP" which tolerates both chiral and achiral forms.
      - A thioester linkage, given by the pattern "[#6](=O)S", is found and its sulfur is embedded 
        in the above CoA connector.
      - The acyl chain attached at the thioester carbonyl (and not part of the CoA scaffold) 
        is composed solely of carbon atoms and is linear (no branches).
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        (bool, str): Tuple (True, reason) if the molecule is classified as 3-substituted propionyl-CoA(4-);
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
    
    # --- 2. Check for the adenine moiety ---
    # SMARTS for adenine (the fused purine ring): n1cnc2cncnc12
    adenine_smarts = "n1cnc2cncnc12"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine):
        return False, "Adenine moiety not found"
    
    # --- 3. Check for the CoA connector (pantetheine unit) ---
    # Use a SMARTS that tolerates chiral centers by using [CX4](O) instead of C(O)
    coa_connector_smarts = "SCCNC(=O)CCNC(=O)[CX4](O)C(C)(C)COP"
    coa_connector = Chem.MolFromSmarts(coa_connector_smarts)
    coa_matches = mol.GetSubstructMatches(coa_connector)
    if not coa_matches:
        return False, "CoA connector (pantetheine fragment) not found"
    # Collect all atom indices from any connector match for cross-checking
    coa_connector_atoms = set()
    for match in coa_matches:
        coa_connector_atoms.update(match)
    
    # --- 4. Check for the thioester linkage and its attachment ---
    # Thioester pattern: a carbonyl carbon bonded to an oxygen (=O) and sulfur (S)
    thioester_smarts = "[#6](=O)S"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    ts_matches = mol.GetSubstructMatches(thioester)
    if not ts_matches:
        return False, "Thioester linkage ([#6](=O)S) not found"
    
    # Ensure that at least one of these thioester links has its sulfur atom in the CoA connector
    valid_ts = []
    for match in ts_matches:
        carbonyl_idx, sulfur_idx = match[0], match[1]
        if sulfur_idx in coa_connector_atoms:
            valid_ts.append(match)
    if not valid_ts:
        return False, "No thioester linkage found whose sulfur is part of the CoA connector"
    
    # --- 5. Validate the acyl chain attached at the thioester carbonyl ---
    # The acyl chain is defined as the part attached to the carbonyl (except the S)
    def extract_linear_chain_length(carbonyl_atom, exclude_idx):
        # Identify carbon neighbors of carbonyl besides the one with index exclude_idx (the S).
        neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != exclude_idx]
        if not neighbors:
            return None  # No acyl chain branch found.
        # Pick the first neighbor as the starting carbon of the acyl chain
        start = neighbors[0]
        length = 1  # count first carbon
        prev = carbonyl_atom
        current = start
        while True:
            # Get the next carbon neighbor excluding where we came from
            nxt = [nbr for nbr in current.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev.GetIdx()]
            if len(nxt) == 0:
                break  # end of chain reached
            elif len(nxt) == 1:
                length += 1
                prev, current = current, nxt[0]
            else:
                # Branching detected
                return None
        return length

    acyl_chain_valid = False
    acyl_chain_length = None
    for ts in valid_ts:
        carbonyl_idx, sulfur_idx = ts[0], ts[1]
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        # Get neighbor carbons of the carbonyl not the sulfur
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors() 
                          if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != sulfur_idx]
        if not acyl_neighbors:
            continue    # no acyl chain attached
        chain_len = extract_linear_chain_length(carbonyl, sulfur_idx)
        if chain_len is None:
            continue  # branched chain found: invalid
        acyl_chain_valid = True
        acyl_chain_length = chain_len
        break
    if not acyl_chain_valid:
        return False, "Acyl chain is branched or contains non-carbon atoms"
    
    # Optionally, impose a minimum acyl chain length (set threshold as desired; here threshold is 3)
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