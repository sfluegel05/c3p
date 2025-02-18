"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
But wait, in the example SMILES, the phosphate is connected via an oxygen to a carbon. So this SMARTS would match that. Then, after confirming the presence of the phosphate ester, we need to check if the carbon attached to the oxygen is allylic (adjacent to a double bond).

So the code would:

1. Check for the presence of the phosphate ester group (C-O-P(=O)(O)O).
2. For each match, get the carbon attached to the oxygen.
3. Check if that carbon is adjacent to a double bond.

For example, in the SMILES ...C(C)=C\COP(O)(O)=O, the carbon attached to the oxygen (O in COP) is connected to another carbon (C in C(C)=C). That carbon has a double bond. So the carbon attached to the oxygen is adjacent to a double bond.

So the code would:

- Find all matches of the phosphate ester pattern.
- For each match, get the carbon atom (the C in [C]-O-P...).
- Check if that carbon has any neighboring atoms that are part of a double bond.

Wait, the carbon attached to the oxygen (let's call it C1) is connected to other atoms. If any of C1's neighbors are part of a double bond, then C1 is allylic.

Alternatively, check if C1 has a neighbor that is part of a double bond. For example, in the structure C-O-P..., the C is connected to O and possibly other carbons. If any of those other carbons are part of a double bond, then C is allylic.

Wait, the allylic position is a carbon adjacent to a double bond. So if the C attached to the oxygen (C-O-P) is adjacent to a double bond, then it's allylic.

But in the example, the structure is C=C-C-O-P. So the C attached to O is the third carbon from the double bond. Wait, no: in ...C(C)=C\COP..., the C attached to O is the one after the double bond. So the structure is C=C-C-O-P. The C attached to O is adjacent to the double bond (the previous carbon is part of the double bond). So the C attached to O is adjacent to a carbon that is part of a double bond. So that C is allylic.

Wait, the allylic position is the carbon adjacent to the double bond. So if the C attached to O is adjacent to a carbon that is part of a double bond, then it's allylic.

So the code should check if the C attached to O has a neighbor (any neighbor) that is part of a double bond.

For example, for the C in C-O-P, check each neighbor of C. If any neighbor is part of a double bond, then it's allylic.

So in code:

matches = mol.GetSubstructMatches(phosphate_ester_pattern)
if not matches:
    return False, "No phosphate ester group"

for match in matches:
    # The first atom in the match is the carbon attached to the oxygen
    c_idx = match[0]
    carbon = mol.GetAtomWithIdx(c_idx)
    # Check all neighbors of this carbon
    for neighbor in carbon.GetNeighbors():
        # Check if the bond between carbon and neighbor is a double bond
        bond = mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx())
        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
            # The carbon is adjacent to a double bond (it's part of a double bond)
            # But wait, this would mean the carbon itself is part of a double bond, which is not the case.
            # We need to check if the neighbor is part of a double bond with another atom.
            for bond2 in neighbor.GetBonds():
                if bond2.GetBondType() == Chem.BondType.DOUBLE and bond2.GetBeginAtomIdx() != c_idx and bond2.GetEndAtomIdx() != c_idx:
                    # The neighbor is part of a double bond with another atom
                    return True, "Phosphate ester attached to allylic carbon"

    # Alternatively, check if any neighbor of the carbon has a double bond
    for neighbor in carbon.GetNeighbors():
        for bond in neighbor.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and (bond.GetBeginAtom() == neighbor or bond.GetEndAtom() == neighbor):
                # The neighbor is part of a double bond
                return True, "Phosphate ester attached to allylic carbon"

return False, "Phosphate not attached to allylic carbon"

This approach checks if any neighbor of the carbon attached to the phosphate's oxygen is part of a double bond. If so, then the carbon is allylic.

Additionally, check for the polyprenol chain's length and number of double bonds. For example, require at least three double bonds in the molecule.

Putting it all together:

from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol_phosphate(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphate ester group (C-O-P(=O)(O)O)
    phosphate_ester_pattern = Chem.MolFromSmarts("[C]-O-P(=O)([O])[O]")
    if not mol.HasSubstructMatch(phosphate_ester_pattern):
        return False, "No phosphate ester group"

    # Find the carbon attached to the phosphate's oxygen
    matches = mol.GetSubstructMatches(phosphate_ester_pattern)
    for match in matches:
        c_idx = match[0]
        carbon = mol.GetAtomWithIdx(c_idx)
        # Check if this carbon is adjacent to a double bond
        for neighbor in carbon.GetNeighbors():
            # Check if the neighbor has any double bonds
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Check if the bond involves the neighbor and another atom
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetIdx() != c_idx:
                        # Found a double bond adjacent to the carbon
                        # Now check the polyprenol chain has multiple double bonds
                        double_bonds = sum(1 for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
                        if double_bonds >= 3:
                            return True, "Phosphate ester attached to allylic carbon of polyprenol chain"
                        else:
                            return False, f"Only {double_bonds} double bonds, need at least 3"

    return False, "Phosphate not attached to allylic carbon"

Wait, but this code has a problem: once it finds a phosphate ester and a double bond adjacent to the carbon, it returns True. But the example might have multiple phosphate groups. Also, the check for double bonds >=3 is inside the loop, which may not be correct. Maybe the double bonds count should be checked outside, after confirming the phosphate is in the correct position.

Revised code:

def is_polyprenol_phosphate(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphate ester group (C-O-P(=O)(O)O)
    phosphate_ester_pattern = Chem.MolFromSmarts("[C]-O-P(=O)([O])[O]")
    if not mol.HasSubstructMatch(phosphate_ester_pattern):
        return False, "No phosphate ester group"

    # Check if the phosphate is attached to an allylic carbon
    allylic_found = False
    matches = mol.GetSubstructMatches(phosphate_ester_pattern)
    for match in matches:
        c_idx = match[0]
        carbon = mol.GetAtomWithIdx(c_idx)
        for neighbor in carbon.GetNeighbors():
            # Check if neighbor has any double bonds
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetIdx() != c_idx:
                        allylic_found = True
                        break
            if allylic_found:
                break
        if allylic_found:
            break

    if not allylic_found:
        return False, "Phosphate not attached to allylic carbon"

    # Check for multiple double bonds in the molecule (indicative of polyprenol chain)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 3:
        return False, f"Only {double_bonds} double bonds, need at least 3 for polyprenol chain"

    # Optional: Check chain length? Not sure, but examples have long chains
    # For example, undecaprenyl has 11 isoprene units (55 carbons?), but maybe not necessary to check here

    return True, "Phosphate ester attached to allylic carbon of polyprenol chain with multiple double bonds"

This code first checks for the presence of a phosphate ester group. Then, it checks if the carbon attached to the phosphate's oxygen is adjacent to a double bond (allylic). Finally, it checks that there are at least three double bonds in the molecule, which is a rough indicator of a polyprenol chain.

Testing this code against the examples:

Take "undecaprenyl dihydrogen phosphate" SMILES: CC(C)=CCC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\COP(O)(O)=O

- The phosphate ester is present.
- The carbon attached to the oxygen (C in COP) is adjacent to a double bond (the C=C in ...C(C)=C\COP...).
- The molecule has multiple double bonds (more than 3).

So this should return True.

Another example: "geranyl phosphate" SMILES: CC(=CCC/C(=C/COP(O)(=O)O)/C)C

- Phosphate ester present.
- The carbon attached to O is part of a C=C (allylic position).
- Double bonds: 2 (maybe), which would fail the >=3 check. But geranyl has two isoprene units (3 double bonds?), wait, geranyl is C10H17OPO3, with two double bonds. So this might return False, but according to the examples, it's a valid polyprenol phosphate. So this code would give a false negative.

Hmm, this is a problem. Geranyl phosphate has two double bonds, but according to the examples, it's a valid polyprenol phosphate. So requiring at least three double bonds is incorrect.

So the check for double bonds >=3 is wrong. The definition says "polyprenol", which is a polymer of isoprene units. The minimum number of isoprene units for a polyprenol is typically three (triprenyl), but sometimes even two (geranyl) is considered. So perhaps the double bond check should be >=2.

Alternatively, remove the double bond count check and rely on the structure. But that might lead to false positives.

Alternatively, the problem is that geranyl phosphate is a prenol phosphate, but not a polyprenol phosphate. Wait, the user's examples include "geranyl phosphate" as a valid example. So according to the problem statement, it should be considered a polyprenol phosphate. But geranyl has two isoprene units (prenol is a single isoprene unit?), so maybe the definition of polyprenol requires at least three. But the user's examples include it, so the code must accept it.

Thus, the double bond check should be >=1. Or perhaps not check the number of double bonds at all, relying solely on the structure of the phosphate attachment.

But that might lead to false positives. For example, a molecule with a single double bond next to a phosphate would be considered a polyprenol phosphate, which might not be correct.

But according to the problem statement, the examples include geranyl phosphate, which has two isoprene units (maybe two double bonds?), so the code must accept it. Hence, the double bond check should be adjusted to >=1.

So modifying the code to check for at least one double bond:

    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 1:
        return False, "No double bonds, not a polyprenol chain"

But even this might not be sufficient. Alternatively, remove the double bond check and just rely on the allylic position and phosphate ester.

Because the definition is about the formal condensation of the terminal allylic hydroxy group of a polyprenol. The key is the allylic position of the hydroxyl (now replaced by phosphate). So even a single isoprene unit (like geranyl) would qualify if it meets the other criteria.

Thus, perhaps the double bond check is unnecessary. The main criteria are:

1. Presence of a phosphate ester group.
2. The ester oxygen is attached to a carbon that is allylic (adjacent to a double bond).

So the code can omit the double bond count check.

Revised code:

def is_polyprenol_phosphate(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphate ester group (C-O-P(=O)(O)O)
    phosphate_ester_pattern = Chem.MolFromSmarts("[C]-O-P(=O)([O])[O]")
    if not mol.HasSubstructMatch(phosphate_ester_pattern):
        return False, "No phosphate ester group"

    # Check if the phosphate is attached to an allylic carbon
    allylic_found = False
    matches = mol.GetSubstructMatches(phosphate_ester_pattern)
    for match in matches:
        c_idx = match[0]
        carbon = mol.GetAtomWithIdx(c_idx)
        for neighbor in carbon.GetNeighbors():
            # Check if neighbor has any double bonds
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetIdx() != c_idx:
                        allylic_found = True
                        break
            if allylic_found:
                break
        if allylic_found:
            break

    if not allylic_found:
        return False, "Phosphate not attached to allylic carbon"

    return True, "Phosphate ester attached to allylic carbon"

This code checks for the phosphate ester and the allylic position. It should work for geranyl phosphate (which has a double bond adjacent to the phosphate's carbon) and other examples.

Testing geranyl phosphate SMILES: CC(=CCC/C(=C/COP(O)(=O)O)/C)C

The phosphate ester is present. The carbon attached to O is part of the structure C=C/COP. The carbon attached to O is connected to a carbon (C) that is part of a double bond (C=O). Wait, no: the SMILES is /C(=C/COP...)/, so the carbon attached to O is part of a double bond. Wait, the SMILES is CC(=CCC/C(=C/COP(O)(=O)O)/C)C. Let's parse this:

The structure is CC(=CCC/C(=C/COP(O)(=O)O)/C)C. The part /C(=C/COP...)/ indicates a double bond between the two C's. So the carbon attached to O is part of a double bond. Wait, no: the part is C(=C/COP...), which means a carbon connected via double bond to another carbon, which is connected to O-P. So the structure is C=C-OP... So the carbon attached to O is adjacent to a double bond (the C=C). Thus, the code would detect the allylic position.

But in this case, the carbon attached to O is directly part of the double bond. Wait, no: the SMILES is C(=C/COP...), which means the carbon is double-bonded to another carbon, which is connected to O-P. So the carbon attached to O is not part of the double bond but is adjacent to it. So the code would detect that the neighbor (the double-bonded carbon) has a double bond, making the carbon attached to O allylic.

Thus, this code should correctly identify geranyl phosphate as a polyprenol phosphate.

Another example: "all-trans-nonaprenyl diphosphate" SMILES: CC(C)=CCC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\COP(O)(=O)OP(O)(O)=O

This has a phosphate ester and the carbon attached to O is adjacent to a double bond. The code should return True.

Thus, this code should handle the examples correctly.
</think>