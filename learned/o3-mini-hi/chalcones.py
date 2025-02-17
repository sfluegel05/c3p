"""
Classifies: CHEBI:23086 chalcones
"""
#!/usr/bin/env python3
"""
Classifies: Chalcones – compounds containing a 1,3-diphenylpropenone core (benzylideneacetophenone)
or the dihydro derivative, that is, molecules that in one part have an exocyclic carbonyl,
with one substituent being directly part of an aromatic ring and the other substituent continuing via a 2-carbon (unsaturated or saturated) chain ending in an aromatic ring.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule belongs to the chalcone (or chalcone derivative) class based on its SMILES.
    Our approach is to (1) find an exocyclic carbonyl (C=O) on a ketone,
    (2) verify that one side of the carbonyl is directly attached to an aromatic ring (a benzene ring),
    and (3) check that the other side leads via a 2-atom linker (either unsaturated, i.e. CH=CH, or saturated, i.e. CH2CH2)
        to a terminal aromatic ring.
        
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a chalcone (or chalcone derivative), False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper: check unsaturated (α,β-unsaturated) branch.
    def check_unsaturated_branch(alpha, parent_idx):
        # Look for a neighbor of alpha (other than the carbonyl carbon) connected by a double bond.
        for nbr in alpha.GetNeighbors():
            if nbr.GetIdx() == parent_idx:
                continue
            bond = mol.GetBondBetweenAtoms(alpha.GetIdx(), nbr.GetIdx())
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # For the candidate beta = nbr, check if one of its other neighbors is aromatic.
                for betaNbr in nbr.GetNeighbors():
                    if betaNbr.GetIdx() == alpha.GetIdx():
                        continue
                    if betaNbr.GetIsAromatic():
                        return True
        return False

    # Helper: check saturated (dihydrochalcone) branch.
    def check_saturated_branch(alpha, parent_idx):
        # Expect a chain of two single-bonded carbons ending in an aromatic fragment.
        for nbr in alpha.GetNeighbors():
            if nbr.GetIdx() == parent_idx:
                continue
            bond1 = mol.GetBondBetweenAtoms(alpha.GetIdx(), nbr.GetIdx())
            if bond1.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Now from this candidate second atom, search for a single-bond neighbor (not alpha)
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() == alpha.GetIdx():
                    continue
                bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                if bond2.GetBondType() != Chem.BondType.SINGLE:
                    continue
                if nbr2.GetIsAromatic():
                    return True
        return False

    # Go through all atoms and search for a suitable carbonyl center.
    for atom in mol.GetAtoms():
        # Look only for carbon atoms
        if atom.GetAtomicNum() != 6:
            continue

        # Check if this carbon is part of a carbonyl group: i.e. it has at least one double-bonded oxygen neighbor.
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    oxy_neighbors.append(nbr)
        if len(oxy_neighbors) < 1:
            continue  # not a carbonyl

        # In chalcones the carbonyl carbon is exocyclic (not part of a ring).
        if atom.IsInRing():
            continue

        # Now get the two carbon neighbors (excluding the oxygen)
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # We expect a ketone so two carbon substituents.
        if len(carbon_neighbors) != 2:
            continue

        # We require that one of the two neighbors be directly aromatic.
        # Try both assignments (one branch will be the direct aromatic side, the other the linker branch)
        for i in range(2):
            direct_side = carbon_neighbors[i]
            linker_side = carbon_neighbors[1-i]
            if not direct_side.GetIsAromatic():
                continue  # direct side must be aromatic

            # Try unsaturated branch first: expect linker_side -> (via double bond) -> next atom -> aromatic.
            unsat_ok = check_unsaturated_branch(linker_side, atom.GetIdx())
            # Alternatively, try saturated branch: expect linker_side -> (via two single bonds) -> aromatic.
            sat_ok = check_saturated_branch(linker_side, atom.GetIdx())
            if unsat_ok or sat_ok:
                if unsat_ok:
                    branch_desc = "α,β-unsaturated ketone variant with terminal aromatic ring via a vinyl linker"
                else:
                    branch_desc = "dihydrochalcone variant with terminal aromatic ring via an ethylene linker"
                reason = f"Contains chalcone core (direct aromatic group attached to carbonyl and linker branch: {branch_desc})"
                return True, reason

    return False, "Chalcone core not found (expected connectivity: aromatic ring - C(=O) - (CH=CH or CH2CH2) - aromatic ring)"

# You can uncomment the test lines below to try the function:
# test_smiles = "O=C(\\C=C\\c1ccccc1)c1ccccc1"  # trans-chalcone
# result, reason = is_chalcones(test_smiles)
# print(result, reason)