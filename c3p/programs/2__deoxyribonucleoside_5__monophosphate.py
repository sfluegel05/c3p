"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies chemical entities of the class 
2'-deoxyribonucleoside 5'-monophosphate.
Definition: A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    The strategy is heuristic. We require that the molecule contains a phosphate group attached
    (via an oxygen) to a 5-membered sugar ring that (1) has exactly one hydroxyl group on the ring carbons 
    (thus missing the 2'-OH) and (2) is connected to a heterocycle (nucleobase).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches the intended features, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens (to permit identification of hydroxyl groups).
    mol = Chem.AddHs(mol)

    ring_info = mol.GetRingInfo()

    # Flag to note if we have found a suitable deoxyribose sugar ring.
    found_deoxy = False

    # Iterate over atoms looking for phosphorus (P) atoms.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "P":
            continue

        # For each P atom, get neighbor oxygens.
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "O":
                continue
            # We want the oxygen that connects P with a carbon (the sugar carbon)
            # Typically the bond connecting the phosphate to the sugar is a single bond.
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            # Now look for a carbon neighbor of this oxygen that is not the phosphorus.
            for o_nbr in nbr.GetNeighbors():
                if o_nbr.GetIdx() == atom.GetIdx():
                    continue
                if o_nbr.GetSymbol() != "C":
                    continue

                sugar_carbon = o_nbr  # candidate sugar carbon
                # Check if this carbon belongs to a 5-membered ring.
                rings = [ring for ring in ring_info.AtomRings() if sugar_carbon.GetIdx() in ring and len(ring) == 5]
                if not rings:
                    continue  # not in a 5-membered ring

                # For each 5-membered ring found, check if it is sugar-like.
                for ring in rings:
                    # In a typical furanose sugar ring there should be four carbons and one oxygen.
                    atom_symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
                    if atom_symbols.count("O") != 1 or atom_symbols.count("C") != 4:
                        continue

                    # Count hydroxyl (-OH) groups on carbons that are in the ring.
                    oh_count = 0
                    for idx in ring:
                        a = mol.GetAtomWithIdx(idx)
                        if a.GetSymbol() != "C":
                            continue
                        # Check each neighbor of this carbon.
                        for nbr2 in a.GetNeighbors():
                            # We consider a substituent to be a hydroxyl if it is an oxygen with at least one hydrogen.
                            if nbr2.GetSymbol() == "O":
                                # Exclude oxygens that are in the ring itself.
                                if nbr2.GetIdx() in ring:
                                    continue
                                # Check if this oxygen has an explicit hydrogen.
                                # (We added hydrogens so they should be explicit.)
                                for subnbr in nbr2.GetNeighbors():
                                    if subnbr.GetSymbol() == "H":
                                        oh_count += 1
                                        break
                    # In a 2'-deoxyribose, the sugar ring should have only one free hydroxyl (on the 3'-position);
                    # whereas a ribose would have an extra OH at the 2'-position.
                    if oh_count != 1:
                        # If more than one OH is found, then it is more like a ribose.
                        continue

                    # Next, as a rough check of the nucleoside nature we require that one of the ring carbons
                    # is attached (outside the sugar ring) to a heterocycle (assumed to contain at least one nitrogen).
                    has_nucleobase = False
                    for idx in ring:
                        atom_in_ring = mol.GetAtomWithIdx(idx)
                        # Look at neighbors outside the ring.
                        for ext in atom_in_ring.GetNeighbors():
                            if ext.GetIdx() in ring:
                                continue
                            # If the neighbor is a carbon or nitrogen and is in a (non‐aliphatic) ring,
                            # assume it might be part of the nucleobase.
                            if ext.GetSymbol() in ["N", "C"]:
                                # If the atom is in any ring with size > 5 then assume it is aromatic heterocycle.
                                for r in ring_info.AtomRings():
                                    if ext.GetIdx() in r and len(r) > 5:
                                        has_nucleobase = True
                                        break
                                if has_nucleobase:
                                    break
                        if has_nucleobase:
                            break
                    if not has_nucleobase:
                        # Even if a sugar-phosphate fragment is found, without a nucleobase we do not classify it
                        # as a nucleoside monophosphate.
                        continue

                    # If all checks pass then we have found a 2'-deoxyribose sugar ring attached at the 5'-position.
                    found_deoxy = True
                    break  # no need to check further rings for this candidate

                if found_deoxy:
                    break  # sugar candidate found for this phosphate group
            if found_deoxy:
                break
        if found_deoxy:
            break

    if found_deoxy:
        return True, "Found 2'-deoxyribose ring (with one hydroxyl on the ring) attached at 5'-position to a phosphate and a nucleobase"
    else:
        return False, "Did not find the required 2'-deoxyribonucleoside 5'-monophosphate fragment"

# For testing purposes
if __name__ == "__main__":
    test_smiles = [
        "Nc1nc2n([C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3)c(=O)[nH]c2c(=O)[nH]1",  # 8-oxo-dGMP
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)([O-])=O)O2)c(=O)[nH]c1=O"  # dTMP(-)
    ]
    for s in test_smiles:
        res, reason = is_2__deoxyribonucleoside_5__monophosphate(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")