"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: Biflavonoid

Definition: A flavonoid oligomer obtained by the oxidative coupling of at least two 
units of aryl‐substituted benzopyran rings (or their substituted derivatives), resulting 
in the two ring systems being joined together by a single atom or bond.

This heuristic performs the following checks:
  • Basic sanity: the molecule must have a minimum ring count (≥6) and a molecular weight >400 Da.
  • The presence of at least one flavonoid-like substructure unit is determined by matching any
    of three SMARTS patterns that cover flavone-, flavan-, or flavanone-like cores.
  • Overlapping substructure matches (sharing at least one atom) are clustered as candidate units.
  • Finally, we require that at least two candidate flavonoid units are present and that there is a 
    direct linking bond (a single bond between two aromatic atoms) connecting atoms from different units.
Note: This remains a heuristic classification.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is considered a biflavonoid based on its SMILES string.

    Args:
        smiles (str): FLAVOR SMILES string for input molecule.
    
    Returns:
        bool: True if the molecule is considered a biflavonoid, else False.
        str: Reason for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic sanity checks: count rings and check molecular weight
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings < 6:
        return False, f"Too few rings ({num_rings}) to be a biflavonoid"
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a biflavonoid"

    # Define SMARTS patterns for flavonoid-like units.
    # We use three patterns:
    #  1. Flavone-like core (chromen-4-one)
    #  2. Flavan-like core (saturated at C2 of the pyran ring)
    #  3. Flavanone-like core (phenyl attached to a carbonyl-bearing pyran ring)
    smarts_patterns = [
        "c1ccc2c(c1)oc(=O)c(c2)",   # flavone-like unit
        "c1ccc2c(c1)OC(c2)",         # flavan-like unit
        "c1ccc(C2=CC(=O)OC2)cc1"      # flavanone-like unit
    ]
    
    patterns = []
    for sp in smarts_patterns:
        patt = Chem.MolFromSmarts(sp)
        if patt is None:
            return False, f"Error generating pattern from SMARTS: {sp}"
        patterns.append(patt)
    
    # Collect all substructure matches for the defined patterns.
    # Each match is stored as a frozenset of atom indices.
    matches = []
    for patt in patterns:
        for match in mol.GetSubstructMatches(patt):
            matches.append(frozenset(match))
    
    if not matches:
        return False, "No flavonoid-like substructure found"
    
    # Cluster overlapping matches into candidate flavonoid units.
    # We cluster any matches that share at least one atom.
    candidate_units = []
    for m in matches:
        placed = False
        for group in candidate_units:
            if any(len(m & other) >= 1 for other in group):
                group.append(m)
                placed = True
                break
        if not placed:
            candidate_units.append([m])
            
    num_units = len(candidate_units)
    if num_units < 2:
        return False, f"Found only {num_units} flavonoid unit(s); at least two are needed for a biflavonoid"

    # Check for a linking bond between at least two candidate groups.
    # We require a direct single bond between aromatic atoms from different candidate units.
    connection_found = False
    for i in range(len(candidate_units)):
        atoms_i = set().union(*candidate_units[i])
        for j in range(i+1, len(candidate_units)):
            atoms_j = set().union(*candidate_units[j])
            for a in atoms_i:
                for b in atoms_j:
                    bond = mol.GetBondBetweenAtoms(a, b)
                    if bond is not None:
                        # Check that the bond is a single bond and both atoms are aromatic.
                        if bond.GetBondType() == Chem.BondType.SINGLE:
                            atom_a = mol.GetAtomWithIdx(a)
                            atom_b = mol.GetAtomWithIdx(b)
                            if atom_a.GetIsAromatic() and atom_b.GetIsAromatic():
                                connection_found = True
                                break
                if connection_found:
                    break
            if connection_found:
                break
        if connection_found:
            break

    if not connection_found:
        return False, "No linking bond detected between candidate flavonoid units (or bond not between aromatic atoms)"
    
    return True, f"Found {num_units} flavonoid unit(s) connected; molecular weight = {mol_wt:.1f} Da"

# Simple testing with an example from the provided list.
if __name__ == "__main__":
    # Example: (+)-morelloflavone (a known biflavonoid)
    test_smiles = "Oc1ccc(cc1)[C@@H]1Oc2cc(O)cc(O)c2C(=O)[C@H]1c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1"
    res, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid:", res)
    print("Reason:", reason)