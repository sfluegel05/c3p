"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python
"""
Classifies: 7-hydroxyisoflavones
Definition: A hydroxyisoflavone compound having a hydroxy group at the 7-position.
The function is_7_hydroxyisoflavones takes a SMILES string and:
  1. Attempts to find an isoflavone-like bicyclic core.
  2. Checks that one of the aromatic carbons on the fused A-ring (defined here as the set of 
     carbons from the first ring of the SMARTS) carries a free hydroxy (-OH) substituent.
Notes:
  - A free –OH here means that the oxygen is not double-bonded (as in carbonyl) and appears
    to be directly attached to the aromatic carbon with no additional heavy atoms (i.e. not part
    of a glycoside or sulfate).
  - The scheme is heuristic and may not work perfectly on every decorated isoflavone.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    It does so in two main steps:
      1. Checks for an isoflavone-like core. Here we use a SMARTS pattern for a 2-phenylchromen-4-one core:
             c1ccc2c(c1)occ2=O
      2. Among the atoms participating in the first (A-ring) part of the core, we look for a free -OH
         substituent that is directly attached to an aromatic carbon.
         
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a 7-hydroxyisoflavone, False otherwise.
        str: A reason message for the classification.
    """
    # Parse SMILES to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # To allow for explicit hydrogen checking on –OH, add hydrogens
    mol_with_H = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for an isoflavone (2-phenylchromen-4-one) core.
    # The pattern below looks for a fused bicyclic system:
    #   "c1ccc2c(c1)occ2=O"
    # The first (fused) aromatic ring (the A-ring) is given by atoms corresponding to the part "c1ccc2c(c1)"
    core_smarts = "c1ccc2c(c1)occ2=O"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in SMARTS definition for isoflavone core"
    
    core_matches = mol_with_H.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Isoflavone core not found"
    
    # Define a helper function: determines whether the candidate oxygen (attached to aromatic carbon)
    # is a free hydroxy group.
    def is_free_OH(o_atom):
        # The oxygen should be sp3 (typical for –OH) and not part of a double bond (e.g. carbonyl).
        # In the Mol with H, check that it has at least one hydrogen.
        # Also, if it is bonded to any heavy atoms besides the one aromatic carbon,
        # then it is likely substituted.
        # (We assume that during substructure matching no extra substituents show up on the OH.)
        if o_atom.GetAtomicNum() != 8:
            return False
        # Check bonds from oxygen: count only bonds to non-hydrogen atoms (besides the aromatic carbon)
        heavy_neighbor_count = 0
        has_H = False
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                has_H = True
            else:
                heavy_neighbor_count += 1
        # For a free hydroxyl, we expect exactly one heavy neighbor (the aromatic carbon)
        # and at least one hydrogen.
        return (heavy_neighbor_count == 1) and has_H

    # We now loop over each core match.
    # Note: the SMARTS was defined as "c1ccc2c(c1)occ2=O" and so the match tuple (if it matches)
    # returns 7 atoms in order. We assume that the first ring (the A-ring) is represented
    # by the atoms from positions 0,1,2 and 4 (these come from "c1", "c", "c", and "c(c1)" parts).
    #
    # We then inspect each of these candidate aromatic carbons for an externally attached free -OH.
    candidate_indices = [0, 1, 2, 4]
    found_7OH = False
    for match in core_matches:
        for pos in candidate_indices:
            core_c_idx = match[pos]
            aromatic_c = mol_with_H.GetAtomWithIdx(core_c_idx)
            # Verify the atom is aromatic and is a carbon
            if aromatic_c.GetAtomicNum() != 6 or not aromatic_c.GetIsAromatic():
                continue
            # Check neighbors of this aromatic carbon that are not in the core match.
            for nbr in aromatic_c.GetNeighbors():
                # If the neighbor is not part of the matched core,
                # it may be substituent that gives the –OH.
                if nbr.GetIdx() in match:
                    continue
                # We want a substituent that is O and qualifies as a free hydroxyl.
                if nbr.GetAtomicNum() == 8:
                    # Exclude if the bond between aromatic carbon and oxygen is a double bond.
                    bond = mol_with_H.GetBondBetweenAtoms(aromatic_c.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondTypeAsDouble() == 2.0:
                        continue
                    if is_free_OH(nbr):
                        found_7OH = True
                        break
            if found_7OH:
                break
        if found_7OH:
            break

    if found_7OH:
        return True, "Isoflavone core detected with free hydroxy group on the A-ring (candidate 7-OH)"
    else:
        return False, "Isoflavone core found, but no free hydroxy group was detected on the A-ring at the presumed 7-position"

# Example usage (you may remove or comment these out for production):
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "CC(C)=CCc1c(O)cc2occ(-c3ccc(O)cc3O)c(=O)c2c1O",  # luteone
        "CC(C)=CCc1c(O)cc2occ(-c3ccc(O)cc3)c(=O)c2c1O",   # wighteone
        "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O",                 # 7-hydroxyisoflavone
        "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)cc(O)c2c1=O",      # isowighteone
        # False positive test example (e.g. glycosylated - these may show extra oxygen connectivity):
        "C[C@@H]1Cc2c3O[C@@]4(O)[C@@H](C)O[C@H](C)c5c4c4c6c(cc(=O)c(c(O)c2[C@@H](C)O1)c36)oc1c(O)cc(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)c(c41)c5=O"
    ]
    for s in test_smiles:
        result, reason = is_7_hydroxyisoflavones(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")