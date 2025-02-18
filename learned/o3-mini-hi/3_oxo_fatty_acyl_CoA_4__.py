"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: 3-oxo-fatty acyl-CoA(4-)

Definition:
  An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups 
  of any 3-oxo-fatty acyl-CoA.
  
This improved implementation:
  1. Checks that the overall molecule has a net formal charge of -4.
  2. Searches for a combined motif (via SMARTS) that requires the acyl chain to be directly
     connected (through a thioester bond) to a CoA fragment. The motif is:
           "C(=O)CC(=O)SCCNC(=O)CCNC(=O)"
  3. For each occurrence of this motif the code identifies the acyl‐carbonyl atom (the first C in the motif)
     and “walks” into the fatty acyl chain (i.e. the branch that is not part of the attached CoA fragment).
     It then requires that the fatty acyl part (R–) contains at least 10 carbon atoms and that at least one
     carbon–carbon double bond (i.e. unsaturation) is present.
     
These extra checks are intended to avoid false positives where the motif is present superficially but the acyl part
does not have the expected chain length or unsaturation (as in, for example, short-chain or otherwise non‐fatty acyls).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if the molecule (given as a SMILES string) belongs to the class
    3-oxo-fatty acyl-CoA(4-).

    Checks performed:
      1. The molecule must have a net formal charge of -4.
      2. The molecule must contain a combined motif that ensures its 3-oxo fatty acyl chain is 
         directly connected to a CoA fragment. The SMARTS used is:
             "C(=O)CC(=O)SCCNC(=O)CCNC(=O)"
      3. From an atom in the motif, we attempt to extract the acyl chain “tail” (the fatty acid part)
         that branches off from the acyl carbon (the first C of the motif). It then must meet both:
             a. A minimum number of carbon atoms (>=10) 
             b. Presence of at least one C=C (non‐carbonyl) double bond.
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA(4-), False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Check overall formal charge -- must be -4 for CoA(4-)
    if Chem.GetFormalCharge(mol) != -4:
        return False, f"Molecule has formal charge {Chem.GetFormalCharge(mol)} (expected -4)."

    # Combined SMARTS for the key connection motif.
    combined_pattern = Chem.MolFromSmarts("C(=O)CC(=O)SCCNC(=O)CCNC(=O)")
    if combined_pattern is None:
        return False, "Error constructing combined SMARTS pattern."

    # Look for matches of motif.
    matches = mol.GetSubstructMatches(combined_pattern)
    if not matches:
        return False, "The expected 3-oxo fatty acyl-CoA motif is not found."

    # We now define a helper to “walk” from the acyl carbon into the fatty acyl chain.
    # Starting from the acyl carbon neighbor that is not part of the motif, we perform
    # a breadth-first search (BFS) over atoms that are carbons.
    def get_acyl_chain_info(start_atom_idx, motif_indices):
        """Return (n_carbon, has_double) where:
           - n_carbon is the total number of connected carbon atoms (in the fatty acyl tail)
             reachable from start_atom_idx (excluding atoms in motif_indices).
           - has_double is True if at least one bond (between carbons) in that region is a double bond.
        """
        visited = set()
        queue = [start_atom_idx]
        visited.add(start_atom_idx)
        n_double = 0

        while queue:
            a_idx = queue.pop(0)
            atom = mol.GetAtomWithIdx(a_idx)
            for bond in atom.GetBonds():
                nb = bond.GetOtherAtom(atom)
                nb_idx = nb.GetIdx()
                # Only traverse carbon atoms that are not in the motif.
                if nb.GetAtomicNum() == 6 and nb_idx not in motif_indices:
                    if nb_idx not in visited:
                        visited.add(nb_idx)
                        queue.append(nb_idx)
                    # Check if the bond is a double bond (and not a carbonyl, which would have been
                    # part of the motif if from the acyl carbon).
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        n_double += 1
        return len(visited), (n_double > 0)

    # Check each occurrence of the motif.
    THRESHOLD_CARBONS = 10  # minimum number of carbons expected in the acyl chain part
    for match in matches:
        # In the SMARTS "C(=O)CC(=O)SCCNC(=O)CCNC(=O)",
        # the first atom (index 0) is the acyl carbonyl.
        acyl_carbon_idx = match[0]
        acyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
        # Look for neighbor(s) of the acyl carbon that are NOT in the motif.
        acyl_tail_found = False
        for nb in acyl_carbon.GetNeighbors():
            if nb.GetIdx() not in match and nb.GetAtomicNum() == 6:
                # Found a candidate for the fatty acyl tail.
                chain_length, has_double = get_acyl_chain_info(nb.GetIdx(), set(match))
                if chain_length >= THRESHOLD_CARBONS and has_double:
                    acyl_tail_found = True
                    break
        if acyl_tail_found:
            return True, ("Molecule contains the required 3-oxo fatty acyl chain directly connected to a CoA fragment, "
                          "has the correct net formal charge (-4), and the fatty acyl tail has sufficient chain length "
                          "and unsaturation.")
    return False, ("The molecule meets the basic motif and charge criteria, but its acyl chain tail (fatty acyl portion) "
                   "is either too short or lacks a carbon–carbon double bond.")

# Example usage:
if __name__ == "__main__":
    # A sample positive example:
    smiles_example = ("CCCC\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)"
                      "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)"
                      "OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
    print("Result:", result)
    print("Reason:", reason)