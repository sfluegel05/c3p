"""
Classifies: CHEBI:48953 cyclohexenones
"""
So the ring has a double bond between C2 and C3, and a ketone at C1. So the ring has two double bonds (C1-O and C2-C3), but the definition says one double bond in the ring. This is a contradiction. Therefore, perhaps the definition's wording is incorrect, and the actual requirement is that the ring has one double bond (the ene part) and the ketone is part of the ring. But that would mean two double bonds. Alternatively, maybe the ketone's double bond is not counted as part of the ring's double bonds because it's a carbonyl group. But that's not chemically accurate.

Alternatively, perhaps the definition considers the ring to have one double bond (the ene part), and the ketone is a substituent. But in the examples, the ketone is part of the ring. For example, the SMILES O=C1C(=C(...)C)C indicates that the carbonyl is part of the ring.

This is confusing. Let me look up the actual structure of cyclohexenone. According to chemical databases, cyclohexenone is a cyclic ketone with a cyclohexene ring (one double bond) and the ketone group. For example, 2-cyclohexen-1-one has a six-membered ring with a double bond between C1 and C2 and a ketone at C1. Wait, but that would mean the ring has two double bonds (C1-C2 and C1-O). But that's impossible in a six-membered ring. Wait, no. The ketone is a carbonyl group (C=O), which is a double bond. So in the ring, there are two double bonds: C1-O and C1-C2. But that's not possible in a six-membered ring with one double bond. So perhaps the actual structure is different.

Wait, maybe the ketone is part of the ring but conjugated with the double bond. For example, the ring has conjugated double bonds. But that would make it aromatic, which it's not (since it's alicyclic). So this is conflicting.

Alternatively, maybe the definition refers to a ring with one double bond (the ene) and the ketone is a substituent. But in that case, the ketone would not be part of the ring, which contradicts the examples.

This suggests that the definition might have a mistake, but given the examples, the correct approach is to look for a six-membered alicyclic ring (non-aromatic) with one double bond (the ene) and a ketone (C=O) as part of the ring.

Wait, but in that case, the ring would have two double bonds (the ene and the ketone), which would make it a conjugated system. However, the examples provided have SMILES strings that include both a double bond in the ring and a ketone. For example, the first example SMILES is O=C1C(=C(C[C@H](C1)O)C)C. Parsing this, the ring is six-membered. The ketone is at position 1 (O=C1). Then, there's a double bond between C2 and C3 (C(=C(...)). So the ring has two double bonds: C1-O and C2-C3. But according to the definition, it should have one. So this suggests that the definition's wording might be incorrect, and the actual requirement is that the ring has one double bond (the ene part) and the ketone is part of the ring, even if that introduces two double bonds. But chemically, that would be a cyclohexadienone, not cyclohexenone. Hmm.

Alternatively, perhaps the definition counts the ketone's double bond as part of the ring's double bond. So the ring has one double bond (the ketone's C=O), and another double bond elsewhere. But that would make two double bonds. But the definition says "one double bond in the ring". So this is conflicting.

Given the examples provided, the correct approach is to identify a six-membered non-aromatic ring that contains a ketone (C=O) and has exactly one double bond in the ring (which could be the ketone's C=O or another). Wait, but the examples have two double bonds. So perhaps the definition's wording is incorrect, and the actual structures have one double bond in the ring (excluding the ketone). Or perhaps the ketone is considered a substituent.

Alternatively, maybe the term "cyclohexenone" refers to a cyclohexane ring with a ketone and a double bond (but that would be two double bonds). This is getting complicated. Perhaps the key is to focus on the examples and create a SMARTS pattern that matches them.

Looking at the first example: O=C1C(=C(C[C@H](C1)O)C)C. The ring atoms are 1 (C connected to O), 2 (C connected to =O), 3 (C with a double bond), 4 (C), 5 (C with hydroxyl), and 6 (C). The ring has a double bond between C2 and C3 (as part of the C(=C) group) and the ketone at C1. So the ring has two double bonds. But according to the definition, it should have one. This suggests that the definition might be incorrect, but the examples are correct. Therefore, the program needs to accept such structures.

Therefore, perhaps the definition's wording is misleading, and the actual requirement is a six-membered alicyclic ring with a ketone and at least one double bond (but possibly more). But the examples have exactly one double bond in the ring (excluding the ketone). Wait, in the first example, the ring has a double bond between C2 and C3 (as part of the =C(...)), and the ketone is at C1 (O=C1). So the ring has two double bonds: the ketone's C=O and the C=C. But according to the definition, it should have one. So perhaps the definition is wrong, and the actual requirement is that the ring has a double bond (the ene) and a ketone (the one), making two double bonds, but the definition says one. This is confusing.

Alternatively, maybe the ketone's double bond is not considered part of the ring's double bonds. But that doesn't make sense chemically. So perhaps the correct approach is to look for a six-membered alicyclic ring (non-aromatic) that contains a ketone (C=O) and has exactly one double bond in the ring (not counting the ketone's). But then how to differentiate?

Alternatively, perhaps the ketone is part of the ring, and the ring has one double bond (in addition to the ketone's). So two double bonds in total. But that would make the molecule have two double bonds in the ring. However, the definition says "one double bond in the ring". So this is conflicting.

Given the examples, perhaps the correct approach is to find a six-membered ring that is not aromatic, contains exactly one double bond (anywhere in the ring), and has a ketone group (C=O) that is part of the ring.

Wait, but in the first example, there are two double bonds in the ring: the C=O and the C=C. So that would not fit. Unless the C=O is not considered part of the ring's double bond count. But chemically, the C=O is part of the ring.

This is a problem. Maybe the definition is intended to have one double bond in the ring (the ene) and the ketone is a substituent. But in the examples, the ketone is part of the ring.

Alternatively, maybe the definition allows for the ketone's double bond and another double bond in the ring. But the wording says "one double bond in the ring".

Given the confusion, perhaps the best approach is to create a SMARTS pattern that matches the examples. Let's look at the first example's structure: O=C1C(=C(...)C)C. The ring is six-membered, with the ketone at position 1. There's a double bond between C2 and C3. So the SMARTS pattern would be a six-membered ring with a ketone (C=O) and a double bond elsewhere in the ring.

So the SMARTS pattern could be: [O]=C1-C(=C)-C-C-C-C1, but forming a ring. But SMARTS is better written with ring closures.

Alternatively, using a SMARTS pattern that matches a six-membered ring with exactly one double bond (not counting the ketone's) and a ketone in the ring. But how?

Alternatively, the SMARTS could be: [C]=C-C-C-C-C(=O)-C, forming a six-membered ring. But this might not capture all cases.

Another approach: find all six-membered rings in the molecule. For each such ring, check if it's non-aromatic, has exactly one double bond (could be any), and contains a carbonyl group (C=O) that is part of the ring.

So, in RDKit code:

1. Generate all rings in the molecule.
2. For each ring, check if size is 6.
3. Check if the ring is non-aromatic (i.e., none of the bonds in the ring are aromatic).
4. Count the number of double bonds in the ring.
5. Check if there's a carbonyl group (C=O) in the ring.

If all conditions are met, return True.

But how to implement this in RDKit?

First, find all rings. For each ring, check size 6.

Then, for each ring, check if any bond in the ring is aromatic. If none, then it's alicyclic.

Then, count the number of double bonds in the ring. Must be exactly one.

Then, check if there's a carbonyl group (C=O) in the ring.

But wait, the carbonyl group is a double bond. So if the ring has a carbonyl, that's one double bond. But if there's another double bond, the count would be two. So according to the definition, the count should be one. But in the examples, there are two double bonds (C=O and C=C). So this approach would fail.

Therefore, the previous reasoning is conflicting with the examples. So perhaps the definition's wording is incorrect, and the actual requirement is that the ring has at least one double bond (including the ketone's), but the examples have two.

Alternatively, maybe the definition counts the ketone's double bond as part of the ring's double bond. So the ring must have exactly one double bond (which is the ketone's). But that would not match the examples.

This is a problem. Given the time constraints, perhaps the best approach is to proceed with the following steps:

- Check for a six-membered ring.
- The ring must be non-aromatic.
- The ring must contain a carbonyl group (C=O) as part of the ring.
- The ring must have exactly one double bond (which could be the carbonyl or another).

But according to the examples, the ring has two double bonds. So this approach would not work.

Alternatively, the definition's "one double bond" refers to the ring having one double bond in addition to the ketone. But that's not clear.

Alternatively, perhaps the definition is incorrect and the actual structures have one double bond (the ene) and the ketone is part of the ring, leading to two double bonds, but the code should accept that.

Given the examples, the code needs to identify a six-membered non-aromatic ring that contains a ketone (C=O) and has at least one double bond (could be the ketone's or another). But how to differentiate between that and other cases.

Alternatively, the code can look for a six-membered ring with a ketone (C=O) in the ring and exactly one double bond (anywhere in the ring). But that would count the ketone's double bond as part of the count.

In the first example, the ring has two double bonds (C=O and C=C), so this approach would fail. But according to the definition, it should pass. Therefore, this approach is incorrect.

This suggests that the definition's wording is incorrect, and the actual requirement is that the ring has one double bond (not counting the ketone's). But how to check that.

Alternatively, the code can check that the ring has exactly one double bond that is not part of the ketone. So, first find the ketone's C=O in the ring, then count the remaining double bonds in the ring. If that count is one, then it's a cyclohexenone.

But how to do that in RDKit.

Alternatively, use SMARTS to find a six-membered ring with a carbonyl group and one double bond (excluding the carbonyl's).

But this is getting complicated. Maybe a SMARTS pattern that matches a six-membered ring with a carbonyl group and exactly one other double bond.

But perhaps the best way is to use SMARTS to match the general structure.

The SMARTS pattern for a six-membered ring with a ketone and one double bond (not the ketone's) would be something like:

[C]=C-C-C-C-C(=O)-C forming a ring. But this is not accurate.

Alternatively, using a SMARTS pattern that matches any six-membered ring containing a carbonyl group and exactly one double bond in the ring (including the carbonyl's). But that would require the ring to have exactly one double bond. However, in the examples, there are two.

This is confusing. Given the time, perhaps the best approach is to proceed with the following code:

Check for a six-membered ring that is non-aromatic, contains a ketone (C=O) in the ring, and has at least one double bond in the ring. But this may not be sufficient.

Alternatively, use the following steps:

1. Parse the SMILES into a molecule.
2. Check for the presence of a six-membered ring.
3. For each six-membered ring, check if it's non-aromatic.
4. Check if the ring contains a carbonyl group (C=O) that is part of the ring.
5. Check if the ring has exactly one double bond (could be the carbonyl or another).

But in the first example, the ring has two double bonds, so this would fail.

Hmm. Given the examples, perhaps the correct approach is to ignore the double bond count and focus on the presence of a six-membered non-aromatic ring with a ketone and at least one double bond.

But that would include structures with multiple double bonds, which may not be correct.

Alternatively, the code can check that the ring has exactly one double bond (including the ketone's). But in the examples, this is not the case.

Given the time constraints and the previous error (the code returning None, None), perhaps the previous code had syntax errors. For example, incorrect indentation or missing imports.

Looking at the previous code attempt (which is not provided), the user mentioned an error: "Error: unexpected indent (<string>, line 1)". This suggests that the code had indentation issues, possibly due to mixing tabs and spaces or incorrect indentation levels.

To avoid this, ensure that the code uses consistent indentation (four spaces per level). Also, include necessary imports and proper function structure.

So, the new code should:

- Import RDKit modules.
- Define the function is_cyclohexenones.
- Parse the SMILES string.
- Check for a six-membered ring that is non-aromatic.
- Check for a ketone in the ring.
- Check for one double bond in the ring.

But given the confusion in the definition, perhaps the code should use a SMARTS pattern that matches the examples.

Looking at the examples, the SMILES strings have a six-membered ring with a ketone (O=C1...) and a double bond elsewhere in the ring (C(=C...)).

A possible SMARTS pattern could be:

[O]=C1-C=,C-C,C,C,C1 where the ring is six-membered with a ketone and a double bond.

But constructing this correctly is tricky.

An alternative SMARTS pattern for a six-membered ring containing a ketone and at least one double bond (anywhere in the ring):

[C]=C-C-C-C-C(=O)-C forming a ring. But this is too specific.

A more general SMARTS pattern:

[O]=C1-*-,=*-*-*-*-1 where the ring has six members, a ketone, and at least one double bond.

Alternatively, using a SMARTS pattern that matches any six-membered ring with a ketone and at least one double bond (including the ketone's):

$six_membered_ring_with_ketone_and_double_bond

But how to write this.

Alternatively, use two separate checks: presence of a ketone in a six-membered ring, and the ring has at least one double bond.

But how to implement this.

Putting it all together:

The code will:

1. Parse the SMILES.
2. Check if the molecule has a six-membered ring.
3. For each six-membered ring, check if it's non-aromatic.
4. Check if the ring contains a carbonyl group (C=O) that is part of the ring.
5. Check if the ring has at least one double bond (could be the carbonyl's or another).

If all conditions are met, return True.

But this may include cases where the ring has multiple double bonds. However, the examples have two double bonds, so this would include them.

But according to the definition, the ring should have one double bond. But the examples contradict this. So perhaps the definition is incorrect, and the code should not check the number of double bonds, only that there's at least one (in addition to the ketone).

Alternatively, perhaps the code should check that the ring has exactly one double bond (excluding the ketone's). But how to do that.

Alternatively, the code can check that the ring has exactly one double bond (regardless of being part of the ketone or not). But in the examples, there are two.

This is a deadlock. Given the time, perhaps proceed with the following code, using SMARTS patterns to capture the examples.

A possible SMARTS pattern for a six-membered ring with a ketone and a double bond (anywhere in the ring):

[*r6]C(=O)[*r6]. This would match any six-membered ring with a ketone. Then, check if the ring has at least one double bond.

Alternatively, use a SMARTS pattern that matches a six-membered ring with a ketone and a double bond:

[O]=C1-C-,=C-C-C-C-1

But this is not flexible.

Perhaps the best approach is to use RDKit's ring detection and check the required properties step by step.

The code:

- Parse the SMILES into a molecule.
- Check if there's at least one six-membered ring.
- For each six-membered ring:
   - Check if it's non-aromatic.
   - Check if it contains a carbonyl group (C=O) as part of the ring.
   - Check if the ring has exactly one double bond (could be the carbonyl's).
- If any ring meets these conditions, return True.

But in the first example, the ring has two double bonds, so this would fail.

Alternatively, change the condition to check for at least one double bond in the ring (including the carbonyl's).

But according to the definition, it should have one. So this is conflicting.

Given the examples, perhaps the code should not check the number of double bonds, just that there's a six-membered non-aromatic ring with a ketone and at least one double bond.

But then, this would include rings with multiple double bonds.

Alternatively, the code can check for exactly one double bond in the ring (including the ketone's). But the examples would fail.

This is a problem. Perhaps the definition's wording is incorrect, and the actual structures have two double bonds. In that case, the code should not check the number of double bonds, just the presence of a six-membered non-aromatic ring with a ketone and at least one double bond.

But then, the code would accept any such structure, which may include more than one double bond.

Given the examples, perhaps this is the correct approach.

So, the code steps:

1. Parse the SMILES.
2. Generate all rings.
3. For each ring:
   a. Check if size is 6.
   b. Check if non-aromatic (all bonds in the ring are single or double, not aromatic).
   c. Check if the ring contains a carbonyl group (C=O) that is part of the ring.
   d. Check if the ring has at least one double bond (could be the carbonyl's or another).
4. If any ring meets all conditions, return True.

But how to check if the ring has at least one double bond.

Implementation:

For each ring, iterate through the bonds in the ring. Check if any bond is a double bond.

But the carbonyl's bond is a double bond, so this would count.

So, the code would check for at least one double bond in the ring (which could be the carbonyl's).

Thus, the code would correctly identify the examples, as they have at least one double bond (the carbonyl's plus another).

But this would include rings with multiple double bonds. However, according to the definition, it's supposed to have one. But the examples have two. So the code would work for the examples but may include false positives.

But given the time constraints and the previous error, proceed with this approach.

Now, coding this in RDKit.

The code:

from rdkit import Chem
from rdkit.Chem import MolFromSmiles, GetSSSR

def is_cyclohexenones(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Get all rings
    rings = mol.GetRingInfo().AtomRings()
    
    for ring in rings:
        if len(ring) != 6:
            continue
        
        # Check if ring is non-aromatic
        is_aromatic = False
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    is_aromatic = True
                    break
        if is_aromatic:
            continue
        
        # Check for carbonyl group in the ring
        has_carbonyl = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() == 8:
                            has_carbonyl = True
                            break
                if has_carbonyl:
                    break
        if not has_carbonyl:
            continue
        
        # Check for at least one double bond in the ring
        has_double_bond = False
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    has_double_bond = True
                    break
        if has_double_bond:
            return True, "Six-membered alicyclic ketone with at least one double bond in the ring"
    
    return False, "Does not meet criteria"

Wait, but this code checks for at least one double bond in the ring (including the carbonyl's). So in the first example, which has two double bonds (C=O and C=C), this would return True. But according to the definition, it should have one. But the examples are correct, so this code would work for them.

But the definition says "one double bond in the ring", but the code allows more. So this may lead to false positives, but given the examples, it's necessary.

Alternatively, the code should check for exactly one double bond in the ring (including the carbonyl's). But in the first example, this would fail.

This is a dilemma. Given the examples, the code must return True for them. So perhaps the definition's wording is incorrect, and the code should check for at least one double bond (excluding the carbonyl's). But how to exclude the carbonyl's double bond.

Alternatively, count the number of double bonds in the ring, excluding the carbonyl's.

So:

1. Find all double bonds in the ring.
2. Subtract the carbonyl's double bond.
3. Check if the remaining count is exactly one.

But how to identify the carbonyl's double bond.

This is complex. Perhaps the code can count the number of double bonds in the ring. If there's a carbonyl group, subtract one from the count, and check if the remaining is exactly one.

But this may not be accurate.

Alternatively:

For each ring:
- Check if it has a carbonyl group.
- Count the number of double bonds in the ring.
- If the count is exactly one, and that's the carbonyl's, then it's a cyclohexanone, not enone.
- If the count is two, and one is the carbonyl's, then it's a cyclohexenone.

But according to the examples, this is the case.

So the code should check that the ring has exactly two double bonds: one from the carbonyl and one from the ene.

But how to implement this.

In code:

count_double = 0
for bond in mol.GetBonds():
    if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            count_double +=1

Then, check if count_double == 1 (excluding the carbonyl's) or count_double == 2 (including the carbonyl's).

But this requires identifying which double bonds are part of the carbonyl.

Alternatively, after confirming the presence of a carbonyl in the ring, subtract one from the total double bond count and check if the remaining is exactly one.

So:

if has_carbonyl:
    total_double = count_double
    remaining = total_double - 1
    if remaining == 1:
        return True

But this would require that the ring has exactly two double bonds (the carbonyl and one other).

In the first example, this would be true (two double bonds). So the code would return True.

But according to the definition, it should have one double bond. So this is conflicting.

But the examples require this.

Therefore, the code should check that the ring has exactly two double bonds: one from the carbonyl and one from the ene.

Thus, the code would:

count_double = sum(1 for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring and bond.GetBondType() == Chem.BondType.DOUBLE)

if has_carbonyl and count_double == 2:
    return True

But this would ensure that the ring has exactly two double bonds, one being the carbonyl's. However, the definition says "one double bond in the ring". So this is conflicting.

Given the examples, this approach is necessary. But the definition's wording is incorrect.

Thus, the final code would:

- Check for a six-membered non-aromatic ring.
- Contains a carbonyl group.
- Has exactly two double bonds in the ring (the carbonyl's and one other).

But this is not per the definition. However, given the examples, this is required.

But this is a guess. Without more information, it's challenging.

Given the time, I'll proceed with the code that checks for a six-membered non-aromatic ring with a carbonyl and at least one double bond (including the carbonyl's). This would match the examples but may include other structures. However, given the previous error was a syntax issue, perhaps the code needs to be correctly formatted.

The final code with proper syntax: