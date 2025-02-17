"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
#!/usr/bin/env python
"""
Classifies: 3-substituted propionyl-CoA(4-)

Definition:
 'An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups
  of any 3-substituted propionyl-CoA; major species at pH 7.3.'

This implementation checks:
  1. That the molecule has a net formal charge of –4.
  2. That the distinctive CoA scaffold is present. In particular we demand that the
     thioester sulfur is linked via the “CoA connector” fragment:
         SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP 
     is found.
  3. That an adenine substructure ("n1cnc2cncnc12") is found.
  4. That a thioester linkage [#6](=O)S is present and that the acyl chain (the part attached
     to the carbonyl on the acyl side) is composed entirely of carbon atoms (i.e. no extra heteroatoms).
     
If any of these tests fail the molecule is not classified as a 3-substituted propionyl-CoA(4-).
"""

from rdkit import Chem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule matches the class of 3-substituted propionyl-CoA(4-) based on its SMILES string.
    
    Expected features:
      - Net formal charge exactly -4.
      - A thioester linkage (pattern [#6](=O)S) connecting an acyl chain to the CoA.
      - A composite CoA "connector" fragment:
            SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP
        which includes the characteristic pantetheine unit.
      - An adenine moiety (with SMARTS "n1cnc2cncnc12").
      - A fatty acyl chain attached via the thioester carbonyl that must be made up solely of carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple with True if classified as 3-substituted propionyl-CoA(4-), otherwise False,
                     and a reason for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Check net formal charge equals -4 ---
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return False, f"Net formal charge is {net_charge} (expected -4)"
    
    # --- 2. Check for the composite CoA scaffold (the key pantetheine connector) ---
    # We remove strict chiral markers to allow for stereochemical variations.
    coa_scaffold_smarts = "SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP"
    coa_scaffold = Chem.MolFromSmarts(coa_scaffold_smarts)
    if not mol.HasSubstructMatch(coa_scaffold):
        return False, "Composite CoA scaffold (pantetheine connector) not found"
    
    # --- 3. Check for the adenine moiety (CoA headgroup) ---
    adenine_smarts = "n1cnc2cncnc12"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine):
        return False, "Adenine moiety not found"
    
    # --- 4. Check for the thioester linkage AND validate the acyl (fatty acyl) chain ---
    # The thioester should match [#6](=O)S.
    thioester_smarts = "[#6](=O)S"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    ts_matches = mol.GetSubstructMatches(thioester)
    if not ts_matches:
        return False, "Thioester linkage ([#6](=O)S) not found"
    
    # We now look at the acyl chain attached on the acyl side of the thioester.
    # For each thioester match (which returns (carbonyl carbon, sulfur) indices), we:
    #   - Identify the neighbor of the carbonyl that is not the sulfur (this is the beginning of the acyl chain).
    #   - Traverse the connected part to ensure that it consists entirely of carbon atoms.
    acyl_chain_valid = False
    for match in ts_matches:
        carbonyl_idx, sulfur_idx = match[0], match[1]
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify neighbor of the carbonyl that is not sulfur and not the carbonyl oxygen.
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors()
                          if nbr.GetIdx() != sulfur_idx and nbr.GetAtomicNum() == 6]
        if not acyl_neighbors:
            continue  # If no proper acyl neighbor is found, check next match.
        start_atom = acyl_neighbors[0]
        # Traverse the connected acyl chain (using DFS) starting from the identified atom.
        visited = set()
        chain_ok = True
        stack = [start_atom.GetIdx()]
        while stack:
            current_idx = stack.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            # If any atom in the chain is not carbon, then reject this thioester match.
            if atom.GetAtomicNum() != 6:
                chain_ok = False
                break
            for nbr in atom.GetNeighbors():
                # Do not traverse back to the carbonyl as that is part of the thioester group.
                if nbr.GetIdx() == carbonyl_idx:
                    continue
                if nbr.GetIdx() not in visited:
                    stack.append(nbr.GetIdx())
        if chain_ok:
            acyl_chain_valid = True
            break
    if not acyl_chain_valid:
        return False, "Acyl chain contains non-carbon atoms; may be an oxygenated or otherwise modified fatty acid"
    
    # --- All tests passed ---
    return True, "Molecule matches structural criteria for 3-substituted propionyl-CoA(4-)"


# Example usage:
if __name__ == "__main__":
    # Example SMILES for tetracosanoyl-CoA(4-)
    smiles_example = (
        "CCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)"
        "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]"
        "([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    )
    result, reason = is_3_substituted_propionyl_CoA_4__(smiles_example)
    print(result, reason)