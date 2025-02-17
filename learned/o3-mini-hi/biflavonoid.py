"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: Biflavonoid

Definition: A flavonoid oligomer that is obtained by the oxidative coupling of 
at least two units of aryl‐substituted benzopyran rings or its substituted derivatives, 
resulting in the two ring systems being joined together by a single atom or bond.

Heuristic improvements:
  - Check for a moderate ring count (≥6) and reasonable molecular weight (>400 Da).
  - Look for two kinds of flavonoid-like substructures with two variants (tight and loose):
       * Tight flavone-like: "c1ccc2c(c1)OC(=O)c(c2)"
       * Loose flavone-like: "c1ccc2c(c1)O*c(=O)*c2"
       * Tight flavan-like:  "c1ccc2c(c1)OC(c2)"
       * Loose flavan-like:  "c1ccc2c(c1)O*c2"
  - Group overlapping matches (if they share >1 atom) into candidate flavonoid units.
  - Require that at least two candidate units exist and that a linking bond or a very small overlap (≤2 atoms) is found
    between them.
Note: This remains a heuristic classification.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES of the molecule.
        
    Returns:
        bool: True if molecule is considered a biflavonoid, else False.
        str: A textual explanation of the outcome.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanity checks: ring count and molecular weight
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings < 6:
        return False, f"Too few rings ({num_rings}) to be a biflavonoid"
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a biflavonoid"
    
    # Define two sets of SMARTS patterns for flavonoid-like units.
    # We use both a 'tight' version and a 'loose' (relaxed) version.
    pattern_smarts = [
        "c1ccc2c(c1)OC(=O)c(c2)",    # tight flavone-like
        "c1ccc2c(c1)O*c(=O)*c2",      # loose flavone-like (allowing a wildcard in the linkage)
        "c1ccc2c(c1)OC(c2)",          # tight flavan-like
        "c1ccc2c(c1)O*c2"            # loose flavan-like
    ]
    
    patterns = []
    for smarts in pattern_smarts:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            return False, f"Error generating pattern from SMARTS: {smarts}"
        patterns.append(pat)
    
    # Collect all substructure matches from all patterns.
    # Each match is stored as a frozenset of atom indices.
    matches = []
    for pat in patterns:
        for match in mol.GetSubstructMatches(pat):
            # Use frozenset so that later we can easily group overlapping matches.
            matches.append(frozenset(match))
    
    if not matches:
        return False, "No flavonoid-like substructure found"
    
    # Group (cluster) overlapping matches into candidate flavonoid units.
    groups = []
    for match in matches:
        added = False
        for group in groups:
            # If match overlaps with any match in the group by more than one atom, group them together.
            if any(len(match & other) > 1 for other in group):
                group.append(match)
                added = True
                break
        if not added:
            groups.append([match])
    num_units = len(groups)
    if num_units < 2:
        return False, f"Found only {num_units} flavonoid unit(s); at least two are needed for a biflavonoid"
    
    # Check that at least one connection exists between the candidate groups.
    connection_found = False
    for i in range(len(groups)):
        group_i_atoms = set.union(*groups[i])
        for j in range(i+1, len(groups)):
            group_j_atoms = set.union(*groups[j])
            common = group_i_atoms & group_j_atoms
            # If they share a very small overlap (up to 2 atoms) treat it as linkage
            if common and len(common) <= 2:
                connection_found = True
                break
            # Otherwise, look for a direct bond between any atoms in the two groups.
            for a in group_i_atoms:
                for b in group_j_atoms:
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
        return False, "Did not detect a linking atom or bond between candidate flavonoid units"
    
    return True, f"Found {num_units} flavonoid unit(s) joined by a short linker; molecular weight = {mol_wt:.1f} Da"

# For basic testing, run one of the provided examples if executed as main.
if __name__ == "__main__":
    # Example: (+)-morelloflavone SMILES
    test_smiles = "Oc1ccc(cc1)[C@@H]1Oc2cc(O)cc(O)c2C(=O)[C@H]1c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1"
    res, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid:", res)
    print("Reason:", reason)