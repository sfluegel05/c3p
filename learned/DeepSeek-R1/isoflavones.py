"""
Classifies: CHEBI:38757 isoflavones
"""
</think>

Okay, let's figure out why the previous code didn't work. The main issue was that all the test cases were missing the core structure. The core_pattern was "c1ccc2c(c1)occc(=O)c2-[c]". Let me visualize that.

The chromen-4-one part is a benzene ring fused to a pyrone (which is a cyclic ketone with an oxygen). The pyrone would have the structure O=C1C=CC2=C1C=CC=C2. The 3-position of the pyrone is where the aryl group attaches. Wait, in the SMARTS pattern, maybe the connectivity is wrong. Let me check the SMARTS notation again.

The original pattern "c1ccc2c(c1)occc(=O)c2-[c]" might not correctly represent the fused benzene-pyrone system. The pyrone part (occc(=O)) might not be correctly structured. Maybe the order of atoms in the ring is incorrect. For example, the ketone oxygen should be part of the pyrone ring. Let me adjust the SMARTS to better represent the 4H-chromen-4-one structure.

Also, the previous code checked for an aryl group connected via a single bond ("-[c]"), but in some isoflavones, the aryl might be directly fused or connected via a double bond? No, the definition says 3-aryl, which is a substituent. So the aryl is attached at the 3-position of the chromenone.

Maybe the SMARTS pattern should be more precise. Let's build the chromen-4-one core correctly. The pyrone ring is O=C1C=C... So the SMARTS for the core should be "[O]=C1C=Cc2c1cccc2" but fused with the benzene. Wait, maybe the correct SMARTS for the chromen-4-one is "c1ccc2c(c1)oc(=O)cc2". Let me draw this:

The benzene ring (c1ccc...) is fused to a 4-membered ring? Wait, no. Chromen-4-one is a benzene fused to a gamma-pyrone. The pyrone is a 4-pyranone, which is a six-membered ring with one oxygen and a ketone. Wait, no. The pyrone is a cyclic ketone with oxygen. The structure is O=C1C=CC=C1. Then fused to a benzene ring. So the correct SMARTS for the chromen-4-one would be "c1ccc2c(c1)oc(=O)ccc2". Let me check that. The benzene (c1ccc...) is connected to the oxygen in the pyrone. The pyrone ring would be O=C1C=CC=C1, but fused to the benzene. So the fused system would be c1ccc2c(c1)oc(=O)ccc2. Wait, maybe the SMARTS should be "c12c(cccc1)oc(=O)ccc2". Let me test this.

Alternatively, using RDKit's SMARTS notation, the core structure is a benzene ring fused to a pyrone where the ketone is at position 4. The correct SMARTS might be "[O]=C1C=Cc2c1cccc2". Wait, perhaps a better approach is to look for the chromen-4-one structure. Let me find a correct SMARTS for that.

Another approach: The chromen-4-one skeleton has the following structure: a benzene ring (c1ccc...) fused to a pyrone ring (O=C1C=CC=C1). The fusion is such that the benzene's two adjacent carbons are part of the pyrone. So the correct SMARTS would be "c1ccc2c(c1)oc(=O)cc2". Let's break this down:

- c1ccc... : benzene ring part.
- c2 connected to c1, then o (oxygen), c(=O) (ketone), then two carbons (cc), and back to c2 to form the pyrone ring.

So the SMARTS pattern for the chromen-4-one core is "c1ccc2c(c1)oc(=O)cc2". Then, the aryl group is attached to the carbon adjacent to the ketone oxygen. Wait, the 3-position in the pyrone ring. Let's count the atoms in the pyrone ring. The pyrone ring would be O=C1-C-C-C-C-1. The numbering starts at the ketone oxygen. So position 3 would be the third carbon in the pyrone ring. In the SMARTS pattern, after the O=C, the next carbon is position 2, then position 3. So the aryl group should be attached to position 3.

In the SMARTS "c1ccc2c(c1)oc(=O)cc2", the pyrone ring is o-c(=O)-c-c. So the third carbon in the pyrone (after O=C) would be the third atom in that ring. So the aryl group would be attached to that carbon. So the SMARTS should have a connection point there.

So modifying the SMARTS to have a substituent at position 3: "c1ccc2c(c1)oc(=O)c([c])cc2". Wait, maybe using a wildcard for the aryl group. Alternatively, the core pattern should be the chromen-4-one with a substituent (aryl) at position 3.

So the correct SMARTS would be "[c]c3oc(=O)c2c1ccccc1c(c2)c3". Hmm, perhaps not. Let me think again. The correct core is the chromen-4-one (benzene fused to pyrone) with an aryl group attached to the third carbon of the pyrone ring.

So the SMARTS pattern should match the chromen-4-one and check that there's an aryl group attached to position 3. Let's build the SMARTS step by step.

The chromen-4-one can be represented as:

c1ccc2c(c1)-oc(=O)ccc2

But the numbering might be off. Alternatively, using RDKit's SMARTS notation, perhaps:

"[O]=C1C=Cc2c1cccc2" for the pyrone part, but fused to the benzene.

Alternatively, the correct SMARTS for the core would be "c1c2c(cccc1)oc(=O)cc2", which represents the fused rings. Then, the third position in the pyrone ring (after the O) would be where the aryl group is attached.

Wait, perhaps the correct approach is to match the chromen-4-one structure and then check that there's a substituent (aryl) attached to the correct position.

Another way: the core is the benzene ring fused to the pyrone, and the aryl group is connected to the carbon adjacent to the oxygen in the pyrone. Let's create the SMARTS for the core structure and the substituent.

The core structure SMARTS could be:

"c1ccc2c(c1)oc(=O)cc2"

Then, the aryl group is attached to the carbon in the pyrone ring that is adjacent to the oxygen. Let's check the atoms in the pyrone ring. The pyrone ring atoms in order would be O, C(=O), C, C, C, C (forming the six-membered ring). Wait, no. The pyrone is a six-membered ring with one oxygen and a ketone. The correct arrangement is O=C1C=CC=C1, but fused to a benzene ring.

Alternatively, the pyrone ring in chromen-4-one is O=C1C=Cc2c1cccc2 (fused to benzene). The position where the aryl is attached (3-aryl) would be the third carbon in the pyrone ring. So in the SMARTS, after O=C1, the next two carbons are positions 2 and 3. The aryl group is attached to position 3.

So the SMARTS for the core would be "c1ccc2c(c1)oc(=O)c([c])cc2". The [c] represents the aryl group attached to position 3.

But SMARTS uses [c] for aromatic carbons. So the pattern would look for an aromatic carbon attached to the pyrone's third position. However, the aryl group could be a substituted benzene or other aromatic ring. So the SMARTS should check that the substituent is part of an aromatic ring.

So the core_pattern should be something like:

"c1ccc2c(c1)oc(=O)c(-[c,c@&H1])cc2"

But maybe more accurately, the substituent should be an aromatic ring. Alternatively, the core_pattern could be the chromen-4-one with any substituent at position 3, and then check if that substituent is part of an aromatic ring.

So the revised code would first check for the chromen-4-one core, then check that the substituent at position 3 is part of an aromatic ring.

In the previous code, the core_pattern was "c1ccc2c(c1)occc(=O)c2-[c]", which might have the oxygen and ketone in the wrong order. The correct order should be "oc(=O)" instead of "occc(=O)".

So fixing the core_pattern to "c1ccc2c(c1)oc(=O)ccc2-[c]" might help. Wait, let's compare:

Original core_pattern: "c1ccc2c(c1)occc(=O)c2-[c]"

This would match a pyrone ring with O followed by three carbons, then the ketone. That's incorrect because the ketone is part of the oxygen's carbon. The correct structure is O=C, so the SMARTS should have "oc(=O)".

So the correct core_pattern should be "c1ccc2c(c1)oc(=O)cc2-[c]". Let's break this down:

- c1ccc2c(c1): benzene ring fused to the pyrone.
- oc(=O): oxygen connected to a carbonyl group.
- cc2: two more carbons to complete the pyrone ring.
- -[c]: the substituent (aryl) attached to the third carbon of the pyrone.

Wait, in the pyrone ring, the atoms would be O, C(=O), then three more carbons? No, the pyrone is a six-membered ring. Let me draw it:

The pyrone ring is O-C(=O)-C-C-C-C, forming a six-membered ring. The benzene is fused to two adjacent carbons in the pyrone. For example, the benzene's two adjacent carbons are part of the pyrone's C-C-C-C part.

Alternatively, the correct SMARTS for the chromen-4-one core with the aryl at position 3 would be:

"c1ccc2c(c1)oc(=O)c(-[c])cc2"

Here, the aryl group (-[c]) is attached to the carbon adjacent to the carbonyl group in the pyrone ring. This would correctly represent the 3-aryl substitution.

So modifying the core_pattern to "c1ccc2c(c1)oc(=O)c(-[c])cc2" might fix the issue.

In the previous code, the core_pattern was "c1ccc2c(c1)occc(=O)c2-[c]", which had the carbonyl in the wrong position. Fixing that to "oc(=O)" instead of "occc(=O)" would correct the core structure.

Once the core is correctly identified, the next step is to verify that the substituent is an aryl group (aromatic ring). The previous code checked if the substituent is part of a 6-membered aromatic ring, but maybe some of the test cases have smaller or larger rings. However, the definition says "3-aryl", which implies an aromatic ring, typically benzene but possibly substituted.

So the revised code should first match the correct core structure, then check that the substituent is part of an aromatic ring (any size? Or at least 6-membered? The examples include some with benzodioxol substituents, like 7-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-1-benzopyran-4-one, which has a benzodioxol group. Benzodioxol is a 5-membered ring with two oxygens. So the aryl group in the definition might include such substituents. Wait, the definition says "3-aryl-4H-chromen-4-one", where aryl is any aromatic substituent. So the substituent can be any aromatic group, not necessarily a benzene ring.

Therefore, the code should check that the substituent is part of an aromatic ring, regardless of size. So the previous code's check for a 6-membered ring was incorrect. That's why many examples were rejected.

In the previous code, after finding the core, it checks if the substituent is part of a 6-membered aromatic ring. But in the test cases, some substituents are benzodioxol (5-membered), which would fail. So the code should instead check if the substituent is part of any aromatic ring, regardless of size.

So the revised code should:

1. Correctly identify the chromen-4-one core with substituent at position 3.
2. Check that the substituent is part of any aromatic ring (could be benzene, benzodioxol, etc.).

So modifying the core_pattern and relaxing the ring checks.

Putting this all together:

The core_pattern should be "c1ccc2c(c1)oc(=O)cc2-[c]". Wait, perhaps better to use a SMARTS that ensures the substituent is at position 3. Let's use:

core_pattern = Chem.MolFromSmarts("[c]1ccc2c(c1)oc(=O)cc2")

Wait, no. The substituent should be attached to the pyrone's third carbon. Let me think again. The correct SMARTS would have the aryl group attached to the carbon adjacent to the carbonyl in the pyrone. So:

The pyrone ring is O=C1-C-C-C-C-C1, fused to benzene. The substituent is at position 3 (third carbon from the carbonyl). So the SMARTS would be:

"c1ccc2c(c1)oc(=O)c(-[c])cc2"

But SMARTS notation requires explicit connections. Alternatively, using a more precise SMARTS:

"c1ccc2c(c1)oc(=O)c([c,c@H1])cc2"

This pattern matches the chromen-4-one core with a substituent (either aromatic or non-aromatic, but likely aromatic) at the correct position.

Once the core is matched, then check if the substituent is part of any aromatic ring.

So in the code, after finding the core matches, iterate over each match and check if the substituent's atom is in any aromatic ring.

Also, the previous code's core_pattern might not have matched because the substituent wasn't properly specified. For example, in formononetin (SMILES: COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O), the core is present. Let's parse its structure. The core is the chromen-4-one (benzene fused to pyrone), with the aryl group (methoxyphenyl) attached at position 3. The SMILES is ...-c1coc2cc(O)ccc2c1=O. The 'coc2' part is the pyrone ring. So the core should be matched by the corrected SMARTS.

In the previous code, the core_pattern didn't match because the SMARTS was incorrect. Fixing the SMARTS to correctly represent the chromen-4-one with substituent at position 3 would allow the examples to be matched.

So the revised code steps:

1. Parse SMILES.
2. Check for the presence of the chromen-4-one core with a substituent at position 3.
3. For each match, check if the substituent is part of any aromatic ring (regardless of size).
4. Return True if any such substituent is found.

Additionally, the previous code's loop over matches might have issues. For example, if there are multiple matches, it should return True if any of them have a valid aryl group.

So the revised code:

- Use the corrected core_pattern.
- Check if the substituent atom is part of any aromatic ring.

Another thing to consider: the substituent could be a single aromatic atom (like a phenyl group) or part of a larger aromatic system. The code should check that the substituent is part of an aromatic ring.

So, in the code:

After finding the core matches, for each match, the substituent's connecting atom is checked. If that atom is in any aromatic ring, then it's considered valid.

Let me adjust the code accordingly.

The core_pattern in SMARTS would be:

"c1ccc2c(c1)oc(=O)cc2-[c]"

But to ensure the substituent is attached to the correct position, perhaps:

"c1ccc2c(c1)oc(=O)c(-[c])cc2"

This would match the chromen-4-one with a substituent (aryl) attached to the carbon adjacent to the carbonyl.

Once this is matched, then check if the substituent's atom (the one attached to the core) is part of an aromatic ring.

In the code:

matches = mol.GetSubstructMatches(core_pattern)
for match in matches:
    # The substituent's connecting atom is the last in the match
    aryl_connector = match[-1]
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if aryl_connector in ring:
            # Check if the ring is aromatic
            is_aromatic = all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
            if is_aromatic:
                return True, "Contains 3-aryl-4H-chromen-4-one skeleton"

Thus, even if the ring is part of a benzodioxol (5-membered), as long as it's aromatic, it's considered valid.

This should fix the previous issues where substituents with non-6-membered aromatic rings were rejected.

Testing this with formononetin's SMILES: COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O.

The core would be matched by the corrected SMARTS. The substituent is the methoxyphenyl group (-c1ccc(OC)cc1), which is aromatic. So the code should return True.

Another example: the SMILES for 7-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-1-benzopyran-4-one. The substituent is a benzodioxol group, which is a 5-membered aromatic ring. The code should detect that the substituent is part of an aromatic ring and return True.

Thus, the revised code should correctly identify the core and the aryl substituent, leading to correct classifications.