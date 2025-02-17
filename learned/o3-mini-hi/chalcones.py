"""
Classifies: CHEBI:23086 chalcones
"""
#!/usr/bin/env python3
"""
Classifies: Chalcones – compounds containing a 1,3-diphenylpropenone core (benzylideneacetophenone)
or the dihydro derivative. In these molecules one part has an exocyclic ketone (C=O) with one substituent
being directly aromatic and the other substituent continuing via a two-atom linker (unsaturated, CH=CH, or
saturated, CH2CH2) to a terminal aromatic ring.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule belongs to the chalcone (or chalcone derivative) class based on its SMILES.
    
    The routine first searches for a carbonyl (C=O) center that is not in a ring and that has exactly two
    carbon neighbors. It then requires that one of these neighbors is aromatic (the “direct” aryl group) and
    that the other neighbor (the “linker side”) is not aromatic but leads via a two-atom chain (either via a double
    bond or via two single bonds) to an aromatic atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a chalcone (or chalcone derivative), False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Helper: check unsaturated (vinyl) branch.
    # Given a candidate branch atom (alpha carbon) attached to the ketone,
    # look for at least one double-bonded neighbor (beta carbon) that is aromatic.
    def check_unsaturated_branch(branch_atom, parent_idx):
        for nbr in branch_atom.GetNeighbors():
            if nbr.GetIdx() == parent_idx:
                continue
            bond = mol.GetBondBetweenAtoms(branch_atom.GetIdx(), nbr.GetIdx())
            # Look for a double bond from alpha (branch_atom) to beta
            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                if nbr.GetIsAromatic():
                    return True
        return False

    # Helper: check saturated (ethylene) branch.
    # Given the branch atom (alpha carbon, non-aromatic) attached to the ketone,
    # look for a single-bonded neighbor (beta carbon) that is aromatic.
    def check_saturated_branch(branch_atom, parent_idx):
        for nbr in branch_atom.GetNeighbors():
            if nbr.GetIdx() == parent_idx:
                continue
            bond = mol.GetBondBetweenAtoms(branch_atom.GetIdx(), nbr.GetIdx())
            if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                if nbr.GetIsAromatic():
                    return True
        return False

    # Iterate over all atoms searching for an exocyclic ketone carbon.
    for atom in mol.GetAtoms():
        # Consider only carbon atoms (atomic number 6)
        if atom.GetAtomicNum() != 6:
            continue

        # Must have at least one double-bonded oxygen (C=O)
        oxy_double = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    oxy_double.append(nbr)
        if len(oxy_double) < 1:
            continue  # not a carbonyl

        # In chalcones the ketone carbon should be exocyclic (not inside a ring)
        if atom.IsInRing():
            continue

        # Get carbon neighbors (exclude oxygen)
        c_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # For a simple ketone we expect exactly two carbon atoms attached.
        if len(c_neighbors) != 2:
            continue

        # Now decide which branch is the “direct aromatic” side and which is the linker side.
        # We require that one substituent is aromatic and the other is not (so that it can be the linking chain).
        for i in range(2):
            direct_side = c_neighbors[i]
            linker_side = c_neighbors[1-i]
            # The direct side must be aromatic.
            if not direct_side.GetIsAromatic():
                continue
            # The linker side must NOT be aromatic (it must serve as the linking atom not part of an aromatic ring).
            if linker_side.GetIsAromatic():
                continue

            # Check linker side for unsaturated (vinyl) possibility
            unsat_ok = check_unsaturated_branch(linker_side, atom.GetIdx())
            # Check linker side for saturated (ethylene) possibility
            sat_ok = check_saturated_branch(linker_side, atom.GetIdx())
            if unsat_ok or sat_ok:
                if unsat_ok:
                    variant = "α,β-unsaturated ketone variant with terminal aromatic ring via a vinyl linker"
                else:
                    variant = "dihydrochalcone variant with terminal aromatic ring via an ethylene linker"
                reason = f"Contains chalcone core (direct aromatic group attached to carbonyl and linker branch: {variant})"
                return True, reason
    return False, ("Chalcone core not found – expected a non‐ring ketone (C=O) with one aromatic substituent "
                   "and one non‐aromatic linker leading (via CH=CH or CH2CH2) to an aromatic ring.")

# Uncomment the lines below to perform a test:
# test_smiles = "O=C(\\C=C\\c1ccccc1)c1ccccc1"  # trans-chalcone
# result, reason = is_chalcones(test_smiles)
# print(result, reason)