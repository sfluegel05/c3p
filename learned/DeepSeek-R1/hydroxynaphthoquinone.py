"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
</think>

Okay, let's see why the previous code didn't work. The main issues were false negatives where valid examples were missed. The error messages say "No naphthalene moiety found" for molecules that should have it. Let me check the naphthalene pattern used.

The previous code used the SMARTS pattern "c1ccc2ccccc2c1", which is a specific arrangement of two fused benzene rings. But naphthoquinones have two ketone groups, so the actual structure might have double bonds in different positions. The rigid pattern might not match when the quinone groups are present, altering the aromatic system.

Ah, right! Naphthoquinones have two ketone groups, which break the aromaticity in those positions. The original naphthalene pattern expects fully aromatic rings, which won't match if the quinone oxygens are present. So the code failed to detect the naphthalene backbone because it's looking for aromatic carbons where there are actually ketones.

To fix this, I need to adjust the naphthalene pattern to account for the quinone groups. Instead of requiring all carbons to be aromatic, the pattern should allow for two ketone groups (C=O) in the naphthalene system. Maybe using a more flexible pattern that includes the quinone structure.

Also, checking for ketones by counting double bonds to oxygen might not be reliable. Using SMARTS to find the quinone groups directly would be better. Let's create a pattern that matches naphthoquinone, which is a naphthalene with two ketone groups. Then, check that this structure exists in the molecule.

So the steps should be:

1. Find a pattern that matches naphthoquinone (two ketones in the naphthalene system).
2. Ensure there's at least one hydroxyl group attached to the naphthoquinone part.

Let me adjust the naphthalene pattern. The naphthoquinone can have different positions for the ketones (like 1,2; 1,4; etc.), so the pattern should account for that. Maybe use a generic naphthalene structure with two ketone groups, regardless of their positions.

Using SMARTS for naphthoquinone: something like C(=O) in two positions of a naphthalene skeleton. But how to represent the fused rings with two ketones?

Alternative approach: Use a SMARTS pattern that matches two fused rings (10 atoms in total) with two ketone groups. For example, something like:

[c]1[c][c][c][c]2[c][c][c][c][c]12 with two of the carbons replaced by C(=O).

But writing that as SMARTS might be tricky. Maybe a better way is to look for two quinone groups in a fused bicyclic system.

Alternatively, find all possible naphthoquinone substructures. For example, a 10-membered bicyclic system with two ketones. Let's create a SMARTS pattern that allows any two adjacent rings with two ketone groups.

Wait, perhaps the correct approach is to first find the naphthoquinone core. A naphthoquinone is a naphthalene derivative where two of the hydrogen atoms are replaced by oxo groups. So the core structure is a bicyclic system (two fused 6-membered rings) with two ketone groups.

So the SMARTS pattern could be:

C(=O)-C1=C-C2=CC=CC=C2C=C1

But that might be too specific. Alternatively, use a more general pattern that allows the two ketones in any position of the fused rings.

Alternatively, split into two parts: check for a naphthalene system (two fused benzene-like rings, possibly with some double bonds) and ensure there are exactly two ketone groups in that system.

Wait, maybe using the original approach but modifying the naphthalene pattern to account for the quinone oxygens. The original code checked for aromatic rings, but in quinones, some carbons are double-bonded to oxygen, breaking aromaticity. So the naphthalene pattern needs to allow for those.

Let me adjust the naphthalene pattern to include possible ketones. For example, the naphthalene skeleton can have two C=O groups. So the SMARTS pattern would be:

c1ccc2c(c1)C(=O)C(=O)cc2

But that's just one possible arrangement. The ketones can be in different positions, like 1,2; 1,4; etc. So perhaps the pattern should be more flexible.

Alternatively, use a SMARTS that allows any two ketone groups in a fused bicyclic system of 10 carbons. Maybe:

C12=C(C(=O)C(=O)C=C1)C=CC=C2

But this is still position-specific. Maybe the correct approach is to look for any bicyclic system (two fused 6-membered rings) with exactly two ketone groups (C=O) in the system.

So, first, find all possible bicyclic systems of 10 atoms (naphthalene-like), then check if there are two ketone groups within that system.

But how to implement that with RDKit? Maybe using ring detection.

Alternative plan:

1. Find all atoms in fused bicyclic systems (two rings, each of 6 members, sharing two adjacent atoms).
2. For each such system, check if there are exactly two ketone groups (C=O) attached to the bicyclic system.
3. Then check if there's at least one hydroxyl group attached to the bicyclic system.

But implementing ring detection for fused bicyclic systems might be complex. Alternatively, use a SMARTS pattern that matches the naphthoquinone structure more accurately.

Wait, perhaps the correct SMARTS for any naphthoquinone would be:

[O]=C=O attached to a fused bicyclic system. Wait, no. Naphthoquinone has two ketone groups (C=O) in the naphthalene system. So the SMARTS pattern should look for a bicyclic system (two fused benzene-like rings) with two C=O groups.

But since the ketones break aromaticity, the pattern needs to account for that. Maybe:

C12=C(C(=O))C(=O)C=CC=C2C=CC=C1

But this is too rigid. Let me look for a more general approach.

Another idea: Use RDKit's ring detection to find atoms in a fused bicyclic system (two rings of 6 atoms each), then check if there are exactly two ketone groups in that system.

But how to code that?

Alternatively, use a SMARTS pattern that matches any two adjacent rings (each of 6 atoms) with two ketone groups. For example:

[O]=C1-C=CC=C2C(=O)C=CC=C12

But this is specific to 1,2-naphthoquinone. Not general enough.

Alternatively, use a SMARTS pattern that allows any two ketone groups in a fused bicyclic system. For example:

[#6]1-,=2-,=3-,=4-,=5-,=6-[#6]2-,=7-,=8-,=9-,=10-,=1[#6])([#8]=[#6])([#8]=[#6])

Hmm, this is getting complicated. Maybe a better approach is to first find all atoms in a fused bicyclic system (two 6-membered rings) and then check if there are exactly two ketone groups in those atoms.

Here's how to approach it:

1. Find all atoms in fused bicyclic systems (two rings of 6 atoms each).
2. For each such system, check if there are exactly two C=O groups in those atoms.
3. Also check that there's at least one hydroxyl group attached to any of those atoms.

But coding this requires complex ring analysis. Let's see if RDKit can help with that.

RDKit's GetRingInfo() can give information about rings. We can iterate over all rings and find pairs that are fused (share at least two atoms). Then, check if the fused rings form a bicyclic system of 10 atoms (since two 6-membered fused rings share two atoms, total 10 atoms).

Once the bicyclic system is identified, check for two ketone groups and at least one hydroxyl attached.

Alternatively, use SMARTS patterns that can capture the naphthoquinone structure with two ketones in any position.

Wait, perhaps the correct SMARTS for naphthoquinone is:

c1c(C(=O))c(C(=O))ccc1

But this is not correct for the fused system. Let me think of a SMARTS pattern that matches a naphthoquinone. For example:

A possible SMARTS for 1,4-naphthoquinone is:

O=C1C=CC(=O)C2=C1C=CC=C2

But this is position-specific. To make it general, perhaps:

[#6]1-,=2-,=3-,=4-,=5-,=6-[#6]2-,=7-,=8-,=9-,=10-,=1[#6])([#8]=[#6])([#8]=[#6]

Not sure. Alternatively, use a SMARTS that allows any two ketone groups in a fused bicyclic system.

Alternatively, use two separate ketone groups anywhere in the molecule, but ensure they are part of a naphthalene-derived system.

Alternatively, the problem with the previous code was that the naphthalene pattern was too strict. The original code looked for a fully aromatic naphthalene, but in naphthoquinones, the ketone groups break the aromaticity in those positions. So the code failed to detect the naphthalene backbone because some carbons are not aromatic.

So the solution is to adjust the naphthalene pattern to allow for the ketone groups. Instead of looking for aromatic carbons, look for a fused bicyclic system with two ketone groups, regardless of aromaticity.

So, the SMARTS pattern for the naphthoquinone core could be:

C12(=O)C=CC=C2C(=O)C=C1

But this is still specific to certain positions. Maybe a better approach is to use a SMARTS pattern that matches any two ketone groups in a fused bicyclic system of 10 atoms.

Alternatively, split the problem into two steps:

1. Check for the presence of two ketone groups (C=O) in the molecule.
2. Check that these ketones are part of a fused bicyclic system (naphthalene-like structure).

But how to check that the ketones are part of the fused system?

Another idea: Use a SMARTS pattern that matches a bicyclic system (two fused rings) with two ketone groups. For example:

[#6]1-,=2-,=3-,=4-,=5-,=6-[#6]2-,=7-,=8-,=9-,=10-,=1[#6])(C(=O)))(C(=O))

But I'm not confident in this pattern. Maybe using RDKit's SMARTS generator isn't the best approach here.

Alternative plan:

- Find all C=O groups in the molecule.
- Check if any two of them are part of a fused bicyclic system (two 6-membered rings sharing two atoms).
- If yes, then check for hydroxyl groups attached to the bicyclic system.

But this requires complex analysis. Let's try to implement this.

First, find all ketone groups (C=O):

ketone = Chem.MolFromSmarts('[CX3]=[OX1]')
ketone_matches = mol.GetSubstructMatches(ketone)
if len(ketone_matches) < 2:
    return False, "Less than two ketone groups"

Then, check if any two ketones are part of a fused bicyclic system.

But how to check if two ketones are in a fused bicyclic system?

Alternatively, for each ketone group, check if the carbon is part of a fused bicyclic system.

But this might be complicated. Let's look for a SMARTS pattern that represents a fused bicyclic system with two ketones.

Alternatively, use the following SMARTS pattern for naphthoquinone:

[C]=[C]-[C]=[C]-[C]=O.[C]=[C]-[C]=[C]-[C]=O

But this is not fused. Hmm.

Wait, perhaps the correct SMARTS for any naphthoquinone is:

C(=O)-C1=C-C(=O)-C2=C-C=CC=C12

But this is for 1,4-naphthoquinone. To generalize, perhaps:

[#6]1(=O)-[#6]=[#6]-[#6]2(=O)-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]12

But this is very specific. Maybe it's better to use a more flexible approach.

Another approach: The original code failed because the naphthalene pattern required aromatic rings, which are not present in the quinone parts. So, perhaps instead of looking for aromatic rings, look for the fused bicyclic structure without considering aromaticity.

The naphthalene pattern in the previous code was "c1ccc2ccccc2c1", which is aromatic. Let's change it to a non-aromatic version, allowing for double bonds and ketones.

The fused bicyclic system can be represented as:

C1=CC=C2C=CC=CC2=C1

But this is a non-aromatic naphthalene. However, in naphthoquinone, there are two ketone groups. So perhaps the SMARTS pattern should be:

C1(=O)=CC=C2C(=O)=CC=CC2=C1

But this is specific to 1,4 positions. To allow any positions, maybe use a more generic pattern.

Alternatively, split the problem into two parts:

1. Check that the molecule has a fused bicyclic system (two 6-membered rings sharing two adjacent atoms).
2. Within that system, there are exactly two ketone groups.
3. At least one hydroxyl group attached to the bicyclic system.

Implementing this requires checking ring systems. Let's try to code this.

First, find all rings in the molecule. Then, find pairs of rings that are fused (share two atoms). For each such pair, check if they form a 10-atom system (since two 6-membered rings sharing two atoms have 6+6-2=10 atoms). Then, check if there are exactly two ketone groups in that system.

Here's how to code this:

mol = Chem.MolFromSmiles(smiles)
ri = mol.GetRingInfo()
rings = ri.AtomRings()

# Find fused ring pairs
fused_pairs = []
for i in range(len(rings)):
    for j in range(i+1, len(rings)):
        intersection = set(rings[i]) & set(rings[j])
        if len(intersection) >= 2:  # Fused rings
            fused_pairs.append( (rings[i], rings[j]) )

# Check each fused pair for being a 10-atom bicyclic system
valid_systems = []
for ring1, ring2 in fused_pairs:
    combined = set(ring1) | set(ring2)
    if len(combined) == 10 and len(ring1) == 6 and len(ring2) ==6:
        valid_systems.append( (ring1, ring2, combined) )

if not valid_systems:
    return False, "No fused bicyclic system found"

# Now check each valid system for exactly two ketones and at least one hydroxyl
for ring1, ring2, combined in valid_systems:
    ketone_count = 0
    for atom_idx in combined:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() !=6:
            continue
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() ==8:
                    ketone_count +=1
    if ketone_count !=2:
        continue  # Check next system

    # Check for hydroxyl groups attached to the combined system
    hydroxyl_found = False
    for atom_idx in combined:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() ==8 and neighbor.GetTotalNumHs() >=1:
                hydroxyl_found = True
                break
        if hydroxyl_found:
            break
    if hydroxyl_found:
        return True, "Fused bicyclic system with two ketones and hydroxyl group"

return False, "No valid hydroxynaphthoquinone system found"

But this code is more complex. Let's adjust the previous code to use this approach.

However, this might be error-prone. Alternatively, use SMARTS patterns that capture the naphthoquinone structure with two ketones and a hydroxyl.

Another idea: Use a SMARTS pattern that matches a bicyclic system with two ketones and at least one hydroxyl group attached. For example:

[OH]-[c]1[c][c][c][c]2[c](=O)[c][c][c][c]12

But this is still specific. Maybe use:

[OH]-[#6]@[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6](=O)~[#6](=O)

Not sure. Alternatively, create a SMARTS pattern that allows any two ketones in a fused bicyclic system with a hydroxyl attached.

Alternatively, use the following steps:

1. Find all ketone groups (C=O).
2. Check if two of them are part of a fused bicyclic system (two 6-membered rings).
3. Check for hydroxyl groups attached to that system.

But implementing this requires checking the rings for each ketone.

Alternatively, the problem in the previous code was that the naphthalene pattern didn't account for the ketone groups. So the solution is to adjust the naphthalene pattern to include the ketones.

For example, the naphthoquinone core can be represented as:

C1=CC=C2C(=O)C=CC(=O)C2=C1

But this is 1,4-naphthoquinone. To make it general, perhaps use a SMARTS pattern that allows any two ketones in the fused system.

A possible SMARTS pattern for any naphthoquinone could be:

[#6]1-,=2-,=3-,=4-,=5-,=6-[#6]2-,=7-,=8-,=9-,=10-,=1[#6].[#8]=[#6].[#8]=[#6]

But this is not correct. Another approach is to use a SMARTS pattern that matches a fused bicyclic system with two ketones:

[#6]1(=O)-[#6]=[#6]-[#6]2(=O)-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]12

But again, this is position-specific.

Given the time constraints, perhaps the best approach is to adjust the naphthalene pattern to allow for ketones by using non-aromatic carbons and checking for the two ketones in the system.

So, modify the original code's naphthalene pattern to be a non-aromatic fused bicyclic system and then check for two ketones and hydroxyl groups.

The previous code's naphthalene pattern was "c1ccc2ccccc2c1", which is aromatic. Let's change it to a non-aromatic version, allowing any connectivity as long as it's a fused bicyclic system.

For example, the SMARTS pattern could be:

C1=CC=C2C=CC=CC2=C1

But this is still a specific arrangement. Alternatively, use a more flexible pattern that matches any two fused 6-membered rings, regardless of aromaticity.

The SMARTS for any two fused 6-membered rings would be:

C1CCCC2C1CCCC2

But this is a decalin structure, which is fully saturated. Not applicable here.

Alternatively, use a SMARTS pattern that matches any fused bicyclic system with 10 atoms:

[#6]@[#6]@[#6]@[#6]@[#6]@[#6]@[#6]@[#6]@[#6]@[#6]

But this is too vague. Alternatively, use the following approach:

1. Find all atoms in a fused bicyclic system (two 6-membered rings).
2. Check if there are exactly two ketone groups in those atoms.
3. Check for hydroxyl groups attached to those atoms.

Implementing this in code:

mol = Chem.MolFromSmiles(smiles)
ri = mol.GetRingInfo()
rings = ri.AtomRings()

# Find all fused pairs of 6-membered rings
valid_systems = []
for i in range(len(rings)):
    ring1 = rings[i]
    if len(ring1) !=6:
        continue
    for j in range(i+1, len(rings)):
        ring2 = rings[j]
        if len(ring2) !=6:
            continue
        intersection = set(ring1) & set(ring2)
        if len(intersection) >=2:
            combined = set(ring1) | set(ring2)
            valid_systems.append(combined)

if not valid_systems:
    return False, "No fused bicyclic system found"

# Check each valid system for ketones and hydroxyl
for system in valid_systems:
    ketone_count =0
    for atom_idx in system:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() !=6:
            continue
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() ==8:
                    ketone_count +=1
    if ketone_count !=2:
        continue

    hydroxyl_found = False
    for atom_idx in system:
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() ==8 and neighbor.GetTotalNumHs() >=1:
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break
    if hydroxyl_found:
        return True, "Fused bicyclic system with two ketones and hydroxyl group"

return False, "No valid hydroxynaphthoquinone found"

This code checks for any fused bicyclic system (two 6-membered rings) with exactly two ketone groups and at least one hydroxyl attached. This should cover the naphthoquinone structure regardless of the positions of the ketones and hydroxyl groups.

Now, integrating this into the function:

The previous code failed because it used an aromatic naphthalene pattern which didn't match the non-aromatic parts due to ketones. By checking for fused rings and counting ketones within them, this approach should be more accurate.

Testing with the example SMILES:

For example, 2,8-dihydroxy-1,4-naphthoquinone (SMILES: C1=CC=C(C2=C1C(C=C(C2=O)O)=O)O) has a fused bicyclic system with two ketones and hydroxyl groups. The new code should detect the fused system, count two ketones, and find the hydroxyl groups.

Thus, the revised code should address the false negatives by correctly identifying the naphthoquinone core even when aromaticity is broken by ketones.

Also, the false positives in the previous attempt were due to incorrect detection of naphthalene systems. By requiring fused bicyclic systems with exactly two ketones and hydroxyl groups, this should reduce false positives.

Finally, adjust the code to use this new approach.