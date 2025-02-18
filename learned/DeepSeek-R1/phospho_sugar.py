"""
Classifies: CHEBI:33447 phospho sugar
"""
</think>

Okay, let's figure out why the previous code didn't work. The main issue seems to be that the sugar ring detection failed for many valid examples. Let's look at the previous approach.

The original code used a SMARTS pattern for the sugar ring: "[O;R]@;=,#[C;R]@;=,#[C;R]@;=,#[C;R]@;=,#[C;R]@;=,#[O;R]". Hmm, that looks too specific. Maybe it's requiring two oxygen atoms in the ring, but some sugars might have different arrangements, like furanoses (5-membered) or pyranoses (6-membered) with varying oxygen positions. Also, the original pattern might not account for different ring sizes or substitution patterns.

Another problem: some examples like beta-D-arabinofuranose 5-phosphate have a 5-membered ring (furanose) with one oxygen in the ring. The original pattern required two oxygens in the ring? Let me check. The pattern starts and ends with O;R, so it's looking for a ring with two oxygens. But many sugars, like ribose (in RNA/DNA) have a single oxygen in the ring (the ring oxygen). Oh, right! The typical sugar ring (like pyranose or furanose) has one ring oxygen. So requiring two oxygens in the ring is incorrect. That's probably why all those examples failed—they have a single O in the ring.

So the first fix is to adjust the sugar pattern to look for a 5- or 6-membered ring with at least one oxygen (the ring oxygen) and multiple hydroxyl groups attached to the ring carbons.

Next, the phosphate detection. The original code used "[C!$(*=O)][OX2][P](=[OX1])([OX2-])[OX2]". Wait, maybe the charge handling is off. The phosphate group might be in different protonation states. Also, the connection to the sugar: the oxygen connecting to P should be part of an ester, i.e., connected to a carbon that's not a carbonyl. The previous check for "[CX4][OX2][P]" might be too restrictive. Let me think: the alcoholic oxygen is connected to a carbon (CX4, which is a carbon with four single bonds), then to P. That part might be okay.

Possible steps to fix:

1. Revise the sugar ring detection to find 5- or 6-membered rings with at least one ring oxygen and multiple hydroxyl groups.
2. Use a more flexible SMARTS for the phosphate ester, accounting for different protonation states (like O-P=O with varying O charges).
3. Ensure the phosphate is attached via an ester linkage to the sugar's hydroxyl.

Let's tackle the sugar ring first. A typical monosaccharide has a 5- or 6-membered ring with one ring oxygen. The hydroxyl groups are on the carbons adjacent. So a better SMARTS might look for a ring of size 5 or 6, containing an oxygen, and with multiple -OH groups attached to the ring carbons.

For example, a pattern like:
- Ring of 5 or 6 atoms, at least one of which is O.
- At least two hydroxyl groups (-OH) attached to ring carbons.

But how to express that in SMARTS? Maybe:

For a 5-membered ring (furanose):
"[C;R5][OH]... etc." But perhaps a better approach is to look for a ring (size 5 or 6) with at least one O in the ring, and at least two -OH groups on ring carbons.

The SMARTS could be something like:
"[$([O;R5]),$([O;R6])]" to find rings with O in 5- or 6-membered rings. Then check that there are at least two hydroxyl groups attached to the ring.

Alternatively, use a combined pattern. Let's try a SMARTS that matches a ring (size 5 or 6) with at least one O in the ring, and at least two OH groups on adjacent carbons.

But maybe splitting into two parts: first find the ring, then check for hydroxyls.

Alternatively, a SMARTS pattern that matches a sugar-like ring:

For furanose (5-membered, one O):
"O1[C@H]([OH])[C@H]([OH])[C@H]([OH])[C@H]1" – but this is too specific. Need a more general approach.

Perhaps a better approach is:

1. Find all rings in the molecule.
2. Check if any ring is 5- or 6-membered and contains an oxygen.
3. For those rings, check that there are at least two hydroxyl groups attached to the ring carbons.

But implementing this with SMARTS might be tricky. Alternatively, use RDKit's ring detection and iterate through the rings.

But the code needs to be a single function without loops, perhaps using SMARTS.

Alternative SMARTS for sugar ring: "[C;R5,C;R6][OH]", but that's not sufficient. Maybe a more flexible approach.

Wait, the original problem was that the SMARTS pattern required two oxygens in the ring. Let's adjust that. Let's create a pattern that matches a 5- or 6-membered ring with at least one oxygen in the ring, and multiple hydroxyl groups attached.

So the sugar pattern could be:

For a 5-membered ring with one O and multiple OHs:
"O@[C@@H]1[C@H](O)[C@H](O)O[C@H]1" – but again, too specific.

Alternatively, a generic pattern that matches any 5- or 6-membered ring with at least one O and at least two OH groups attached to the ring.

Maybe using a SMARTS like:

"([C;R5,C;R6]-[OH]){2,}.(O;R5,O;R6)"

Which would match two or more hydroxyls attached to ring carbons (5 or 6 membered) and an oxygen in the ring. But SMARTS doesn't allow counting like that. Alternatively, use a recursive SMARTS.

Alternatively, use two separate substructure matches:

First, check for a 5- or 6-membered ring with at least one O in the ring.

Then, check that there are at least two -OH groups attached to the ring carbons.

But how to do that in SMARTS?

Alternatively, use two separate checks:

1. Presence of a 5- or 6-membered ring containing an oxygen atom.

2. At least two hydroxyl groups (-OH) attached to ring carbons.

But implementing this in RDKit would require more code.

Alternatively, combine these checks into a single SMARTS.

Let me think. For a 5-membered ring with one O and two OHs:

The SMARTS could be something like:

"[C;R5,C;R6]-[OH].([C;R5,C;R6]-[OH]).O;R5,O;R6"

But SMARTS doesn't work like that. Maybe using logical operators.

Alternatively, use a SMARTS that matches a ring (size 5 or 6) with at least one O in the ring and at least two OH groups on ring carbons.

But constructing such a SMARTS is challenging. Perhaps a better approach is to split into two parts:

First, check for any 5- or 6-membered ring that includes an oxygen.

Then, in the molecule, check that there are at least two hydroxyl groups attached to the carbons in such rings.

But how to do that without iterating through atoms?

Alternatively, use a SMARTS that looks for two hydroxyl groups attached to ring carbons in a ring with oxygen.

But maybe this is getting too complicated. Let's look for existing approaches.

Upon reflection, the original code's sugar_pattern was incorrect because it required two oxygens in the ring, which is not the case for most monosaccharides. The correct approach is to look for a ring with one oxygen (the ring oxygen) and several hydroxyl groups.

So, let's adjust the sugar_pattern to find a 5- or 6-membered ring with at least one oxygen in the ring. Then check for hydroxyl groups.

But how to do that in SMARTS?

A possible SMARTS for a furanose (5-membered ring with one O) is "O1C[C@H](O)C[C@H]1O", but that's too specific. We need a more general pattern.

Alternatively, to find any 5- or 6-membered ring containing an oxygen atom, regardless of substituents:

For 5-membered ring with O: "[O;R5]"

For 6-membered ring with O: "[O;R6]"

So the combined pattern would be "[O;R5,O;R6]". Then, check that the molecule has such a ring.

But this would match any ring with O in 5 or 6 membered ring. Then, additionally check that there are at least two hydroxyl groups attached to the ring carbons.

But how to ensure that?

Perhaps:

After confirming the presence of a ring with O, check that there are at least two -OH groups attached to the ring carbons.

But implementing this in code would require:

1. Find all atoms in 5- or 6-membered rings that are oxygen (ring oxygen).
2. For each such ring, check the attached carbons for hydroxyl groups.

But this requires iterating through rings, which may be complex.

Alternative approach: Use a SMARTS that matches a carbon in a 5- or 6-membered ring with an -OH group. Then count at least two such carbons.

So the SMARTS would be "[C;R5,R6]-[OH]". If there are at least two matches, then it's a sugar.

So combining these steps:

First, check if there's a 5- or 6-membered ring containing an oxygen.

Then, check that there are at least two hydroxyl groups attached to carbons in such rings.

But how to implement this in code?

Here's a possible plan:

1. Check if the molecule has any 5- or 6-membered rings containing oxygen.

   - SMARTS for ring oxygen: "O;R5" or "O;R6"

   So create a pattern like "[O;R5,O;R6]".

   If the molecule has a match for this, proceed.

2. Check that there are at least two -OH groups attached to carbons in 5- or 6-membered rings.

   - SMARTS: "[C;R5,R6]-[OH]"

   Count the number of matches. Need at least two.

So in code:

- Check if the molecule has at least one O in a 5- or 6-membered ring.

- Then, check that there are at least two [C;R5,R6]-OH groups.

But how to do this without using loops? Maybe using two separate SMARTS matches.

So revising the code:

First, check for the presence of a ring oxygen in 5- or 6-membered ring.

ring_oxygen = MolFromSmarts("[O;R5,O;R6]")

if not mol.HasSubstructMatch(ring_oxygen):

    return False, "No sugar ring oxygen detected"

Then, check for at least two hydroxyls on ring carbons:

hydroxyl_on_ring = MolFromSmarts("[C;R5,R6]-[OH]")

hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_on_ring)

if len(hydroxyl_matches) < 2:

    return False, "Insufficient hydroxyl groups on sugar ring"

This approach would handle the sugar detection better.

Now, the phosphate ester check. The previous code used "[C!$(*=O)][OX2][P](=[OX1])([OX2-])[OX2]". But phosphate groups can be in different protonation states. For example, the SMILES might have P(=O)(O)(O) or P(=O)([O-])[O-], etc. The original SMARTS might not capture all variants.

A better approach is to look for the ester linkage: a carbon connected to an oxygen, which is connected to a phosphorus with at least one double bond to oxygen.

So the SMARTS could be "[CX4]-[OX2]-P(=O)". This would match a carbon (with four bonds) connected to an oxygen (single bond), connected to a phosphorus with a double bond to oxygen.

But also, the phosphate could have other oxygen groups. So the pattern could be "[CX4]-O-P(=O)(-O)-O". But allowing for different substituents.

Alternatively, use a SMARTS like "[C!$(*=O)]-[OX2]-P(=O)([OX2-])[OX2]" to capture the ester linkage.

Wait, the key is that the phosphate is esterified to an alcoholic oxygen. So the oxygen connected to P must be connected to a carbon that's not part of a carbonyl group.

So the SMARTS for the phosphate ester could be "[CX4]-[OX2]-P(=O)([OX2-,-])[OX2-,-]". This would allow for different charges.

But in RDKit's SMARTS, how to represent possible charges or different oxygen states?

Perhaps the pattern should be "[CX4]-O-P(=O)(O)(O)", but allowing for deprotonated oxygens (with negative charges).

But SMILES often represents phosphate groups as P(=O)(O)(O), but in reality, at physiological pH, they are deprotonated. However, SMILES strings can vary. For example, the example SMILES for ATP has "OP(O)(=O)OP(O)(=O)OP(O)(O)=O".

So perhaps the SMARTS should look for O-P(=O)(O) groups, but allowing for any number of oxygen atoms connected, as long as there's a P=O group.

The key part is the presence of a phosphate group (P=O) connected via an ester oxygen (O connected to a carbon).

So the SMARTS could be "[CX4]-O-P(=O)". This would match any carbon (not part of a carbonyl) connected to an oxygen, which is connected to a phosphorus with a double bond to oxygen.

But this might miss some cases where the phosphate has other substituents. For example, in a diphosphate like ADP, the phosphate is part of a chain. However, the ester is still the first oxygen attached to the sugar.

So the SMARTS "[CX4]-O-P" would match the ester linkage. Then, check that the phosphorus has at least one double bond to oxygen.

So combining these:

phosphate_ester = MolFromSmarts("[CX4]-[OX2]-P(=O)")

But RDKit's SMARTS allows for checking the presence of a double bond. So this pattern matches a carbon connected to an oxygen (single bond), which is connected to a phosphorus with a double bond to oxygen.

This should capture the phosphate ester group.

Now, putting it all together.

The revised steps are:

1. Check that the SMILES is valid.

2. Check for presence of a 5- or 6-membered ring containing an oxygen atom (the ring oxygen).

3. Check that there are at least two hydroxyl groups attached to the ring carbons.

4. Check for a phosphate ester group: a carbon (not carbonyl) connected via oxygen to a phosphorus with a double bond to oxygen.

Additionally, ensure that the phosphate is attached to an alcoholic oxygen (i.e., the oxygen is connected to a carbon that is not part of a carbonyl group).

So the code would be:

def is_phospho_sugar(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for ring oxygen in 5- or 6-membered ring
    ring_oxygen = Chem.MolFromSmarts("[O;R5,O;R6]")
    if not mol.HasSubstructMatch(ring_oxygen):
        return False, "No sugar ring oxygen detected"
    
    # Check for at least two hydroxyl groups on ring carbons
    hydroxyl_on_ring = Chem.MolFromSmarts("[C;R5,R6]-[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_on_ring)
    if len(hydroxyl_matches) < 2:
        return False, f"Only {len(hydroxyl_matches)} hydroxyls on sugar ring"
    
    # Check for phosphate ester: O-P(=O) connected to non-carbonyl carbon
    phosphate_ester = Chem.MolFromSmarts("[CX4]-[OX2]-P(=O)")
    if not mol.HasSubstructMatch(phosphate_ester):
        return False, "No phosphate ester group found"
    
    return True, "Monosaccharide with phosphate ester"

Testing this against the examples:

Take beta-D-arabinofuranose 5-phosphate: SMILES O[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@@H]1O.

The ring is 5-membered (furanose) with one O. The hydroxyl_on_ring pattern would find multiple OHs. The phosphate ester is attached via O-P(=O), so the SMARTS should match.

Another example: ATP's SMILES includes a sugar ring (ribose) with multiple OHs and a phosphate ester chain.

This revised code should detect the sugar ring and phosphate ester.

Potential issues:

- The hydroxyl_on_ring count might include non-ring hydroxyls if the SMARTS is not specific enough. Wait, the SMARTS is [C;R5,R6]-[OH], which means carbons in 5- or 6-membered rings with an OH. So that should correctly count the hydroxyls on the sugar ring.

Another possible issue: some sugars might have substituents other than OH, but the definition requires the phosphate to be on an alcoholic hydroxy group. The current code checks that the phosphate is attached via an alcoholic oxygen (CX4-O-P), which should be correct.

Testing the previous examples that failed:

For 3'-L-arginyl-AMP: the SMILES has a ribose ring (5-membered with O). The hydroxyl_on_ring would have multiple OHs. The phosphate ester is present. So this should be detected.

Similarly, ATP has a ribose ring with OHs and phosphate esters.

Thus, this revised code should address the previous issues by correctly identifying the sugar ring and phosphate ester.