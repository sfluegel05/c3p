"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: Biflavonoid

Definition: A flavonoid oligomer obtained by the oxidative coupling of at least two 
units of aryl‐substituted benzopyran rings (or their substituted derivatives), resulting 
in the two ring systems being joined together by a single atom or bond.

This heuristic checks for:
  • A minimum ring count (≥6) and molecular weight (>400 Da)
  • The presence of at least two flavonoid‐like substructures defined by two SMARTS:
       - Flavone pattern: "c1ccc2c(c1)oc(=O)c(c2)"
       - Flavan pattern:  "c1ccc2c(c1)OC(c2)"
  • Overlapping matches are grouped as one candidate unit and a linking bond between 
    any atoms of different groups is required.
Note: This remains a heuristic classification.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is considered a biflavonoid based on its SMILES string.
    
    Args:
        smiles (str): Input SMILES string.
        
    Returns:
        bool: True if the molecule is considered a biflavonoid, else False.
        str: Reason for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic sanity checks: ring count and molecular weight.
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings < 6:
        return False, f"Too few rings ({num_rings}) to be a biflavonoid"
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a biflavonoid"
    
    # Define SMARTS patterns for flavonoid-like units.
    # Flavone (chromen-4-one) and flavan cores are used here.
    smarts_patterns = [
        "c1ccc2c(c1)oc(=O)c(c2)",  # flavone-like unit
        "c1ccc2c(c1)OC(c2)"        # flavan-like unit
    ]
    
    patterns = []
    for sp in smarts_patterns:
        pat = Chem.MolFromSmarts(sp)
        if pat is None:
            return False, f"Error generating pattern from SMARTS: {sp}"
        patterns.append(pat)
    
    # Collect all substructure matches for the defined patterns.
    # We store each match as a frozenset of atom indices.
    matches = []
    for patt in patterns:
        for match in mol.GetSubstructMatches(patt):
            matches.append(frozenset(match))
    
    if not matches:
        return False, "No flavonoid-like substructure found"
    
    # Cluster overlapping matches into candidate flavonoid units.
    candidate_units = []
    for m in matches:
        placed = False
        for group in candidate_units:
            # If the current match overlaps with any match in the group by 2 or more atoms, merge it.
            if any(len(m & other) >= 2 for other in group):
                group.append(m)
                placed = True
                break
        if not placed:
            candidate_units.append([m])
            
    num_units = len(candidate_units)
    if num_units < 2:
        return False, f"Found only {num_units} flavonoid unit(s); at least two are needed for a biflavonoid"
    
    # Check for a linking bond between at least two candidate groups.
    connection_found = False
    for i in range(len(candidate_units)):
        atoms_i = set().union(*candidate_units[i])
        for j in range(i+1, len(candidate_units)):
            atoms_j = set().union(*candidate_units[j])
            # Look for a direct bond between any atom in atoms_i and any in atoms_j.
            for a in atoms_i:
                for b in atoms_j:
                    if mol.GetBondBetweenAtoms(a, b) is not None:
                        connection_found = True
                        break
                if connection_found:
                    break
            if connection_found:
                break
        if connection_found:
            break
    
    if not connection_found:
        return False, "No linking bond detected between candidate flavonoid units"
    
    return True, f"Found {num_units} flavonoid unit(s) connected; molecular weight = {mol_wt:.1f} Da"

# Simple testing with an example from the provided list.
if __name__ == "__main__":
    # Example: (+)-morelloflavone (a known biflavonoid)
    test_smiles = "Oc1ccc(cc1)[C@@H]1Oc2cc(O)cc(O)c2C(=O)[C@H]1c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1"
    res, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid:", res)
    print("Reason:", reason)