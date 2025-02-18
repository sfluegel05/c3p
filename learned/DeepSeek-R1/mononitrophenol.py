"""
Classifies: CHEBI:39362 mononitrophenol
"""
</think>

Okay, let's tackle this problem. The goal is to write a Python function using RDKit to determine if a given SMILES string represents a mononitrophenol. From the definition, a mononitrophenol is a phenol with exactly one nitro group substituent. 

First, I need to break down the requirements. Phenol means there's a benzene ring with at least one hydroxyl (-OH) group. Then, there should be exactly one nitro group (-NO2) attached to the ring. The position (ortho, meta, para) doesn't matter as per the definition.

So the steps should be:
1. Check if the molecule has a phenolic hydroxyl group. That means an -OH directly attached to an aromatic ring.
2. Ensure there's exactly one nitro group (-NO2) on the aromatic ring.
3. Make sure there are no other nitro groups elsewhere in the molecule.

Wait, but the nitro group could be anywhere on the ring. Also, the molecule might have other substituents as long as there's only one nitro. For example, the examples include compounds with bromo, methoxy, etc., but only one nitro.

So the approach would be:
- Find all aromatic rings in the molecule.
- For each aromatic ring, check if it has a hydroxyl group attached.
- Then, check if that same ring has exactly one nitro group.
- Also, ensure that the total number of nitro groups in the entire molecule is exactly one.

Wait, but the nitro group must be on the phenolic ring. Because the definition says "phenol carrying a single nitro substituent". So the nitro has to be on the same ring as the hydroxyl.

So the steps in code would be:
1. Parse the SMILES into a molecule.
2. Check for validity.
3. Find all aromatic rings.
4. For each aromatic ring, check if it has at least one hydroxyl (-OH) attached.
5. For those rings that have a hydroxyl, check if they have exactly one nitro group (-NO2) attached.
6. Also, ensure that the entire molecule has only one nitro group in total.

Wait, but maybe the molecule could have multiple rings. For example, a biphenyl system where one ring is phenolic and the other has a nitro. But according to the definition, the nitro must be on the phenol ring. So the nitro must be on the same ring as the hydroxyl.

So the correct approach is:
- Find all aromatic rings that have at least one hydroxyl group attached.
- For each such ring, check if it has exactly one nitro group attached.
- Additionally, ensure that the total number of nitro groups in the molecule is exactly one.

But how to implement this in RDKit?

Using SMARTS patterns might help. For the phenolic hydroxyl: a hydroxyl (-O-) attached to an aromatic carbon. The nitro group is [N+](=O)[O-].

So first, check that there's at least one phenolic -OH. Then, check that there's exactly one nitro group on the same ring as the phenolic -OH. And that the total nitro groups in the molecule is one.

Alternatively, maybe find all nitro groups in the molecule. If there's exactly one, then check if it's attached to an aromatic ring that also has a hydroxyl group.

Let me outline the steps:

1. Parse the SMILES into a molecule. If invalid, return False.
2. Find all nitro groups in the molecule. If count !=1, return False.
3. Find the aromatic ring that the nitro group is attached to.
4. Check if that same ring has at least one hydroxyl (-OH) group attached.
5. If all conditions met, return True.

But how to check if the nitro is on the same ring as the hydroxyl?

Alternatively, find all rings that have both a hydroxyl and a nitro group. If any such ring exists, and the total nitro count is 1.

Let me think of SMARTS patterns.

A phenolic hydroxyl can be matched with [OH]-c where the 'c' is aromatic. But RDKit's SMARTS might need explicit aromatic flags.

A possible SMARTS for a phenolic hydroxyl: [OH]-[c] where [c] is an aromatic carbon.

The nitro group can be matched with [N+](=O)[O-].

But how to ensure they are on the same ring?

Perhaps, create a SMARTS pattern that matches a benzene ring with one -OH and one nitro group.

Like:

c1c([OH])ccc([N+](=O)[O-])c1

But this would only match specific positions. Since the positions can vary, maybe a more general approach is needed.

Alternatively, for each aromatic ring in the molecule, check if it has exactly one -OH and exactly one nitro group. But that's not correct because the definition allows other substituents as long as there's one nitro and the hydroxyl.

Wait, the definition says "phenol carrying a single nitro substituent". So the phenol is the parent, and there's exactly one nitro substituent. The phenol can have other substituents as long as there's one nitro.

Wait no: the definition is "phenol carrying a single nitro substituent". So the phenol (which is a benzene ring with an -OH) has exactly one nitro group attached. The other substituents (like bromo, methoxy, etc.) are allowed as long as the nitro count is exactly one.

Wait, looking at the examples provided, like 4-nitrophenol (Oc1ccc(cc1)[N+](=O)[O-]), which has one nitro. Another example is 2-Methoxy-5-nitrophenol (COc1ccc(cc1O)[N+]([O-])=O), which has a methoxy and a nitro. So the presence of other groups is allowed as long as there's exactly one nitro on the phenolic ring.

So the key points are:
- At least one phenolic -OH (attached to aromatic ring).
- Exactly one nitro group on the same aromatic ring as the phenolic -OH.
- The total number of nitro groups in the molecule is one.

Wait, but the nitro could be on a different ring if the molecule has multiple rings. But according to the definition, the nitro must be a substituent on the phenol. So the nitro must be on the same ring as the hydroxyl.

So the steps would be:

1. Check that the molecule has exactly one nitro group.
2. Check that there's at least one phenolic -OH (i.e., -OH attached to an aromatic carbon).
3. Check that the nitro group is attached to the same aromatic ring as at least one of the phenolic -OH groups.

But how to implement this in RDKit?

Approach:

- Find all nitro groups. If count !=1, return False.
- Find all aromatic rings that have at least one -OH attached.
- For each such ring, check if the nitro group is part of that ring's atoms.
- If any ring satisfies this (has an -OH and the nitro is on it), then return True.

But how to get the ring that the nitro is on?

Alternatively, for the nitro group's atom, find the ring it's part of. Then check if that ring has at least one -OH attached.

Wait, the nitro group is a substituent on the ring. So the nitro's nitrogen is bonded to a carbon in the aromatic ring.

So, for the nitro group's atom (the N), find the neighboring carbon (which is part of the ring), then check if that ring has an -OH.

But in RDKit, how to find the ring that a particular atom is in?

Alternatively, iterate over all rings in the molecule. For each ring, check if it contains the nitro's attachment point (the carbon connected to the nitro group) and also has an -OH group.

Let's outline the code steps:

1. Parse the SMILES into a molecule. If invalid, return False.
2. Check the number of nitro groups. If not exactly one, return False.
3. Find all atoms that are part of nitro groups. For the nitro group, get the attached carbon (the one bonded to the N+).
4. Check if that carbon is in an aromatic ring.
5. For that ring, check if there's at least one -OH group attached to any of its atoms.
6. If all conditions met, return True.

But how to implement step 3-5?

Let's break it down:

First, find the nitro group. A nitro group is [N+](=O)[O-]. So in RDKit, we can use a SMARTS pattern to find nitro groups.

nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')
nitro_matches = mol.GetSubstructMatches(nitro_pattern)
if len(nitro_matches) != 1:
    return False

Each nitro_match is a tuple of atom indices. The nitro group's N is the first atom in the match.

Then, for that N atom, find the neighboring atom (the carbon it's attached to). Because in the nitro group, the N is connected to a carbon (assuming the nitro is attached to the ring via a single bond).

So:

nitro_atom = mol.GetAtomWithIdx(nitro_matches[0][0])
neighbors = nitro_atom.GetNeighbors()
# Assuming the nitro is attached to one carbon (the rest are O)
attached_carbon = None
for neighbor in neighbors:
    if neighbor.GetSymbol() == 'C':
        attached_carbon = neighbor
        break

if not attached_carbon:
    return False, "Nitro group not attached to carbon"

Then, check if the attached_carbon is part of an aromatic ring.

But how to check if it's in an aromatic ring? We can check if the attached_carbon is in any aromatic rings.

Alternatively, get the rings that include the attached_carbon and check if any of them is aromatic.

But RDKit's GetRingInfo() can give the rings. Let's see.

ring_info = mol.GetRingInfo()
aromatic_rings = []
for ring in ring_info.AtomRings():
    is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
    if is_aromatic:
        aromatic_rings.append(ring)

Then, check if the attached_carbon's index is in any of the aromatic rings.

But perhaps a better way is to check if the attached_carbon is aromatic and part of a ring.

attached_carbon_idx = attached_carbon.GetIdx()
for ring in aromatic_rings:
    if attached_carbon_idx in ring:
        # This is the aromatic ring the nitro is attached to
        # Now check if this ring has a hydroxyl group (-OH) attached
        # Iterate over atoms in the ring, check if any has an -OH
        has_hydroxyl = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1 and atom.GetDegree() == 1:
                # Check if the oxygen is part of a hydroxyl group attached to the ring
                # The oxygen should be connected to one carbon (the ring) and have one H
                # Also, the bond should be single
                for bond in atom.GetBonds():
                    if bond.GetBeginAtomIdx() == atom_idx or bond.GetEndAtomIdx() == atom_idx:
                        other_atom_idx = bond.GetOtherAtomIdx(atom_idx)
                        other_atom = mol.GetAtomWithIdx(other_atom_idx)
                        if other_atom.GetSymbol() == 'C' and other_atom.GetIdx() in ring:
                            # Check bond type? Maybe not necessary, but could be single
                            if bond.GetBondType() == Chem.BondType.SINGLE:
                                has_hydroxyl = True
                                break
                if has_hydroxyl:
                    break
        if has_hydroxyl:
            return True, "Phenol with one nitro substituent"

return False, "No phenolic ring with nitro group"

Wait, but this approach might miss some cases. For example, the hydroxyl could be attached to a carbon in the ring, but the oxygen is not in the ring. So the hydroxyl is -O-H attached to a ring carbon.

The SMARTS pattern for a phenolic hydroxyl would be [OH]-c, where c is an aromatic carbon. So perhaps using a SMARTS pattern to find all phenolic -OH groups.

phenol_oh_pattern = Chem.MolFromSmarts('[OH]-[c]')
phenol_oh_matches = mol.GetSubstructMatches(phenol_oh_pattern)

If there are no phenol_oh_matches, return False.

Then, for each phenolic -OH, check if its ring has the nitro group.

Alternatively, for each phenolic -OH group, find the ring it's attached to, and check if that ring contains the nitro's attached carbon.

But this could be complicated.

Alternative approach:

1. Check that the molecule has exactly one nitro group.
2. Check that the molecule has at least one phenolic -OH group.
3. Ensure that the nitro group is attached to the same aromatic ring as at least one of the phenolic -OH groups.

So code steps:

mol = Chem.MolFromSmiles(smiles)
if not mol:
    return False, "Invalid SMILES"

# Check for exactly one nitro group
nitro_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[N+](=O)[O-]'))
if len(nitro_matches) != 1:
    return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"

# Check for at least one phenolic -OH
phenol_oh = mol.GetSubstructMatches(Chem.MolFromSmarts('[OH]-[c]'))
if not phenol_oh:
    return False, "No phenolic hydroxyl group"

# Get the attached carbon of the nitro group
nitro_atom_idx = nitro_matches[0][0]
nitro_atom = mol.GetAtomWithIdx(nitro_atom_idx)
attached_carbon = None
for atom in nitro_atom.GetNeighbors():
    if atom.GetSymbol() == 'C':
        attached_carbon = atom
        break
if not attached_carbon:
    return False, "Nitro group not attached to carbon"

# Check if attached_carbon is in an aromatic ring
attached_carbon_ring = None
ring_info = mol.GetRingInfo()
for ring in ring_info.AtomRings():
    if attached_carbon.GetIdx() in ring:
        # Check if the ring is aromatic
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if is_aromatic:
            attached_carbon_ring = ring
            break
if not attached_carbon_ring:
    return False, "Nitro group not attached to aromatic ring"

# Check if any phenolic -OH is on the same ring
for oh_match in phenol_oh:
    oh_oxygen_idx = oh_match[0]
    oh_carbon_idx = oh_match[1]
    # Check if the carbon is in the attached_carbon_ring
    if oh_carbon_idx in attached_carbon_ring:
        return True, "Phenol with one nitro substituent"

return False, "Nitro and hydroxyl not on the same aromatic ring"

This approach should work. Let's test with the example SMILES:

For 4-nitrophenol (Oc1ccc(cc1)[N+](=O)[O-]):
- One nitro group.
- Phenolic -OH present.
- Nitro's attached carbon is in the aromatic ring, and the phenolic -OH's carbon is the same ring. So returns True.

For 2-Methoxy-5-nitrophenol (COc1ccc(cc1O)[N+]([O-])=O):
- One nitro, phenolic -OH present. Both on the same ring. Returns True.

For a molecule with nitro not on the phenolic ring, like if there's another ring, it should return False.

Another example: 3-nitrotyrosine (N[C@@H](Cc1ccc(O)c(c1)[N+]([O-])=O)C(O)=O). The nitro is on the phenolic ring (the benzene with -OH), so returns True.

Potential issues:
- The nitro group's attached carbon must be in an aromatic ring.
- The phenolic -OH must be on the same ring.

Another possible edge case: if the molecule has multiple aromatic rings, but the nitro is on the same ring as a phenolic -OH, and others have -OH but not on that ring. Should still return True.

Thus, this code should correctly identify mononitrophenols.