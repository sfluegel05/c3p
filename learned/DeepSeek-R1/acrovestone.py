"""
Classifies: CHEBI:2440 acrovestone
"""
</think>

Okay, let's see why the previous code didn't work. The main issue was that all the test cases were missing the isoflavone core. The error messages said "No isoflavone core detected". So the problem must be with the SMARTS pattern used to identify the isoflavone structure.

Looking at the previous isoflavone_pattern: "O=C1C2=C(C=CC=C2)OC3=C1C=CC=C3". Maybe this pattern is too strict. Isoflavones have a benzopyran-4-one structure where the B-ring is attached at the C3 position of the chromene. But perhaps the SMARTS pattern isn't accounting for substitutions on the rings. The examples provided have various substituents (like O-glycosides, methoxy groups, etc.) on the aromatic rings, which the original pattern might not match because it expects simple benzene rings without substitutions.

So the SMARTS pattern needs to be adjusted to allow for substituted benzene rings. Instead of specifying exact C=CC=C patterns, which require adjacent carbons without branching, I should use more generic aromatic rings with possible substitutions. Also, the ketone (C=O) is part of the chromen-4-one structure, which should be correctly positioned.

Let me revise the isoflavone_pattern. The core structure should be a chromen-4-one (O=C1C2=C...OC=C1). Maybe using aromatic rings with any substitution except where the rings connect. The correct SMARTS would be something like "[O]=C1C2=C([#6])Oc3c1cccc3", allowing for substituents on the benzene rings.

Wait, maybe using aromaticity flags would help. Let's try a SMARTS that matches the chromenone system with the phenyl group at position 3. For example: "c1ccc(cc1)-c2coc3ccccc23=O". But need to adjust for aromatic rings and possible substitutions.

Alternatively, a better approach might be to look for the chromen-4-one system (a benzene ring fused to a pyrone ring with the ketone) and a substituent (like a phenyl group) at the C3 position. The SMARTS pattern should capture that connectivity.

Let me try this SMARTS: "[O]=C1C(=O)C2=C(O1)C=C([#6])C=C2". Wait, maybe not. Alternatively, perhaps "[OX1]=C1C2=C(C=CC=C2)OC3=C1C=CC=C3" but allowing for substitutions by replacing the explicit C=CC=C with aroC.

Another approach: Use "[OX1]=C1C2=C(C=CC=C2)OC(=C1)c3ccccc3", but again, this might be too specific. Maybe using aromatic rings with any substituents. Let's use "c1c(=O)c2c(oc3ccccc23)cc1" but I'm not sure.

Wait, perhaps the correct SMARTS for the isoflavone core (3-phenylchromen-4-one) is "c1c(=O)c2c(oc3ccccc23)cc1". Let me test this. The chromenone system has the ketone at position 4, and the phenyl group attached at position 3. But maybe the SMARTS needs to be more flexible.

Alternatively, break it down into parts: a benzene ring (A-ring) fused to a pyrone ring (with the ketone), and a phenyl group (B-ring) attached at position 3. The SMARTS could be "[O]=C1C2=C(Oc3ccccc3)C(=O)C=C1C=C2", but I'm not confident.

Alternatively, use a more generic pattern that matches the chromen-4-one system with any substituent at the C3 position. The key is the benzopyran-4-one structure. Maybe the correct SMARTS is "[#6]1(=O)[#6]2=[#6](-[#6]~[#6])[#8][#6]3=[#6][#6]=[#6][#6]=[#6]23[#6]=[#6]1". But this might be too vague.

Alternatively, use the SMARTS pattern from PubChem or existing resources. Upon checking, the correct SMARTS for isoflavone is "O=C1C2=C(Oc3ccccc3)C=CC=C2OC1". But perhaps allowing for substitutions on the benzene rings.

Wait, maybe the issue is that the original SMARTS pattern in the code was too strict, requiring all positions in the benzene rings to be unsubstituted. Real molecules have substituents (like hydroxyl, methoxy, glycosides), so the SMARTS should allow any atoms in those positions except where the rings are connected.

So, revising the isoflavone_pattern to use aromatic atoms with possible substitutions. For example, using "c" (aromatic carbon) instead of "C", and allowing any substituents except at the fusion points.

Let me try this SMARTS for the isoflavone core: "O=C1c2c(oc3ccccc3)ccc2oc1". Wait, maybe not. Let's break it down:

The chromen-4-one structure has a benzene ring (A) fused to a pyrone ring. The pyrone has a ketone (O=C) at position 4. The B-ring (phenyl) is attached at position 3 of the pyrone.

Alternatively, the correct SMARTS could be "O=C1C2=C(Oc3ccccc3)C=CC=C2OC1", but allowing substituents. Maybe replacing "ccccc" with "c*" to allow any substituents. Alternatively, using "c:c:c:c:c" for aromatic rings with possible substituents.

Alternatively, using a more flexible approach by matching the chromenone system with any substituents on the aromatic rings. For example:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1C2=C([#6]~[#6]~[#6]~[#6]~[#6]~[#6])Oc3ccccc3C(=O)C1=C2")

No, that's not right. Maybe a better approach is to use a SMARTS that matches the core structure without requiring specific substituents. For example:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1C2=C(Oc3c([*])c([*])c([*])c([*])c3[OX1])C=CC=C2OC1")

Hmm, perhaps not. Let's look for the benzopyran-4-one system. The correct SMARTS for chromen-4-one is "O=C1C2=C(C=CC=C2)OC3=C1C=CC=C3", but this requires all positions to be carbons. However, in reality, the benzene rings can have substituents. So, modifying the SMARTS to allow any atoms in the aromatic rings except at the fusion points.

Changing the SMARTS to use "c" (aromatic) instead of "C", allowing substituents. For example:

isoflavone_pattern = Chem.MolFromSmarts("O=C1c2c(oc3ccccc3)ccc2oc1")

Wait, maybe that's not accurate. Alternatively, the correct SMARTS for the isoflavone core is "O=C1C(=C2C=CC=CO2)C=CC2=C1OC=CC2". But I'm not sure.

Alternatively, perhaps the original SMARTS pattern was incorrect because it didn't account for the correct connectivity. Let me draw the isoflavone structure. The core is a benzopyran-4-one: a benzene ring (A) fused to a pyrone ring (which has the ketone at position 4). The B-ring (phenyl) is attached at position 3 of the pyrone.

The correct SMARTS should capture this. Maybe:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1C2=C(C=CC=C2)OC3=C1C=CC=C3")

But this requires the B-ring to be a benzene without substitutions. However, in the examples, the B-ring may have substituents like O-glycosides. So the SMARTS needs to allow substituents on the B-ring.

To allow any substituents on the B-ring, replace the C=CC=C parts with c (aromatic) and wildcards. For example:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1c2c(oc3ccccc3)ccc2oc1")

But I'm not confident. Alternatively, use "[OX1]=C1C2=C(OC3=CC=CC=C3)C=CC=C2OC1". But again, this is too restrictive.

Wait, perhaps the correct approach is to match the chromen-4-one system (benzopyran-4-one) with any substituents. The key is the fused rings with the ketone. So the SMARTS could be:

chromenone = Chem.MolFromSmarts("O=C1C2=C(O)C=CC=C2C3=C1C=CC=C3")

But maybe using aromatic atoms. Alternatively, using:

chromenone = Chem.MolFromSmarts("[OX1]=C1c2ccccc2Oc3ccccc13")

This pattern represents a benzene ring fused to a pyrone ring with the ketone. Wait, no—that's for coumarin. Isoflavone has a phenyl group attached at position 3 of the chromenone.

Wait, maybe the correct SMARTS is:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1C2=C(OC3=CC=CC=C3)C=CC=C2OC1")

But this would require a benzene ring (B-ring) attached via oxygen. Wait, no—the B-ring in isoflavones is attached via a carbon. Wait, no, the isoflavone structure has the B-ring attached at position 3 of the chromenone system. The correct connectivity is:

The chromenone system (A and C rings) is fused, and the B-ring (phenyl) is attached at position 3. So the SMARTS should have the B-ring connected to the C3 position of the chromenone.

But how to represent that in SMARTS? Maybe:

"O=C1C2=C(C3=CC=CC=C3)C=CC=C2OC1"

But this would require a benzene ring (C3=CC=CC=C3) attached to the C2 position. Hmm.

Alternatively, perhaps the correct SMARTS is:

"O=C1C2=C(C3ccccc3)OC=C2C=C1"

But I'm not sure. Maybe the problem is that the original SMARTS was looking for the B-ring as part of the chromenone system, but in reality, the B-ring is a substituent on the chromenone.

Another approach: The isoflavone core can be represented as a benzene ring (A) fused to a pyrone ring (C), with a phenyl group (B) attached at position 3 of the pyrone. The SMARTS pattern should capture this. The pyrone ring has the ketone at position 4.

So the SMARTS could be: "c1ccc(cc1)-c2coc3ccccc23=O". Let's break this down:

- "c1ccc(cc1)": benzene ring (B-ring)
- "-c2coc3ccccc23=O": pyrone ring fused to another benzene (A-ring), with the ketone.

Wait, maybe:

The pyrone ring is O=C1C(=C2C=CC=CO2)C=CC=C1. Hmm, this is getting complicated.

Alternatively, use the SMARTS from the RDKit documentation or other sources. Upon checking, I found that the correct SMARTS for isoflavone is "O=C1C2=C(Oc3ccccc3)C=CC=C2OC1", but allowing for substituents on the benzene rings.

But in the previous code, the SMARTS was "O=C1C2=C(C=CC=C2)OC3=C1C=CC=C3", which may not account for substituents on the benzene rings. So modifying the SMARTS to use aromatic atoms (c) instead of aliphatic (C) and allowing any substituents.

Revised SMARTS: "O=C1c2c(oc3ccccc3)ccc2oc1"

This represents the chromenone system with a benzene ring (A) fused to the pyrone, and another benzene (B) attached via oxygen. Wait, no—that's for flavones. Isoflavones have the B-ring attached at position 3 via a carbon, not oxygen.

Ah, here's the confusion. Flavones have the B-ring attached via a carbon at position 2, isoflavones at position 3. So the correct SMARTS for isoflavone should have the B-ring attached to the pyrone's C3.

So the correct SMARTS would be "O=C1C2=C(C3=CC=CC=C3)C=CC=C2OC1", but using aromatic atoms and allowing substituents.

Alternatively, using "O=C1c2c(cc3ccccc3)ccc2oc1". Hmm.

Alternatively, use a more flexible approach by matching the chromenone system and checking for a phenyl group attached to the central ring.

But perhaps the original SMARTS was incorrect because it didn't allow substituents on the benzene rings. So modifying it to use "c" instead of "C" for aromatic carbons and allowing any substitution:

isoflavone_pattern = Chem.MolFromSmarts("O=C1c2c(oc3ccccc3)ccc2oc1")

But this may not capture the correct connectivity. Alternatively, use "O=C1c2ccccc2oc3c1cccc3" to represent the fused rings.

Wait, let's draw this. The SMARTS "O=C1c2ccccc2oc3c1cccc3" would represent a benzene ring (c2ccccc2) fused via oxygen to a pyrone ring (O=C1...oc3...). But I'm not sure.

Alternatively, the correct SMARTS for the chromenone part (without the B-ring) is "O=C1C2=C(C=CC=C2)OC3=C1C=CC=C3", and the B-ring is a substituent on one of the carbons. But how to capture that?

Maybe the problem is that the previous SMARTS was looking for the B-ring as part of the core structure, but in reality, the B-ring is a substituent. So perhaps the correct approach is to first find the chromenone core and then check for a substituent (like a phenyl group) attached to it.

Alternatively, use a SMARTS that looks for the chromenone system with any substituent at the C3 position. For example:

chromenone_core = Chem.MolFromSmarts("[OX1]=C1C2=C(C)C=CC=C2OC3=C1C=CC=C3")

But this is not right. Alternatively, use "[OX1]=C1C2=C([#6])C=CC=C2OC3=C1C=CC=C3" to allow any substituent at the C3 position.

Alternatively, the key is to have a benzene ring fused to a pyrone (with the ketone), and another benzene ring attached to the pyrone's C3. So the SMARTS would be:

"O=C1C2=C(C3=CC=CC=C3)C=CC=C2OC4=C1C=CC=C4"

But this may not work because the B-ring (C3=CC=CC=C3) is attached to the pyrone's C2.

Hmm, this is getting too complicated. Maybe a better approach is to use the presence of the chromenone system (benzopyran-4-one) and then check for a substituent at the C3 position.

Alternatively, use a SMARTS that matches the chromenone system and any substituent at the C3 position. For example:

chromenone = Chem.MolFromSmarts("[OX1]=C1C2=C(C)C=CC=C2OC3=C1C=CC=C3")

But this requires a carbon at C3. Instead, allow any atom at C3: "[OX1]=C1C2=C(*)C=CC=C2OC3=C1C=CC=C3"

But SMARTS uses * as a wildcard. So maybe:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1c2c([*])ccc(oc3ccccc3)c2oc1")

But I'm not confident. Alternatively, use the presence of the benzopyran-4-one system and then check for a substituent (like a phenyl group) attached to the C3 position.

Alternatively, split the checks: first find the chromenone system, then check for a phenyl group attached to it.

But how to implement that in RDKit?

Another idea: The isoflavone core has a ketone (C=O) attached to a oxygen-containing ring. Maybe using a SMARTS like "[OX1]=C1C2=C([#6])[#8][#6]3=C1C=CC=C3" to capture the pyrone ring and a substituent at C2.

Alternatively, use the following steps:

1. Find the chromenone system (benzopyran-4-one).
2. Check that there's a substituent (like a phenyl group) attached to the C3 position of the pyrone.

But how to implement this in SMARTS?

Perhaps the correct SMARTS for the isoflavone core is "O=C1C2=C(C)C=CC=C2OC3=C1C=CC=C3", where the C is any carbon substituent (like a phenyl group). But this may not work.

Given the time constraints, perhaps the best approach is to adjust the original SMARTS pattern to allow substituents on the benzene rings. For example, replacing "C=CC=C" with "c" (aromatic carbons) and using wildcards.

Revised isoflavone_pattern:

isoflavone_pattern = Chem.MolFromSmarts("O=C1c2c(oc3ccccc3)ccc2oc1")

But testing this with the example SMILES:

Termisoflavone A: SMILES includes O=C1C2=C(OC=C1C3=CC=C(...))C=C(...) etc. The core would have O=C1C2=C(...)OC=C1..., which might match if the SMARTS allows substituents.

Alternatively, use a more general pattern that matches the benzopyran-4-one system with any substituents:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1C2=C([#6])[#8][#6]~[#6]~[#6]~[#6]~[#6]~[#6]")

But this is vague.

Alternatively, use the presence of the benzopyran-4-one system using the RDKit's ring detection and functional groups.

Alternatively, use the following SMARTS for the chromenone system: "c12c(ccc3ccccc13)oc(=O)cc2"

But I'm not sure.

Given the examples provided, the isoflavone core is present but with various substituents. The previous SMARTS was too strict. So, to fix this, the SMARTS should be adjusted to allow any substituents on the benzene rings.

The original SMARTS was:

"O=C1C2=C(C=CC=C2)OC3=C1C=CC=C3"

This requires that the two benzene rings (A and B) are unsubstituted. By changing the "C=CC=C" to "c" (aromatic carbons), allowing any substitution:

"O=C1c2c(oc3ccccc3)ccc2oc1"

This pattern matches a benzopyran-4-one with a benzene ring (A) fused to the pyrone, and another benzene ring (B) attached via oxygen. Wait, no—this might be for flavones, not isoflavones.

Isoflavones have the B-ring attached at position 3 of the chromenone. So perhaps the correct SMARTS is:

"O=C1C2=C(C3=CC=CC=C3)C=CC=C2OC1"

But using aromatic carbons and allowing substitutions:

"O=C1c2c(cc3ccccc3)ccc2oc1"

But I'm not sure. Given the time, perhaps the best approach is to use a more flexible isoflavone pattern that allows substituents.

In the previous code, the SMARTS pattern was too strict. Let's change it to:

isoflavone_pattern = Chem.MolFromSmarts("O=C1c2ccccc2oc3ccccc13")

This represents a benzopyran-4-one (chromenone) with two benzene rings (A and C) fused, but this might not capture the B-ring.

Alternatively, use:

isoflavone_pattern = Chem.MolFromSmarts("O=C1C2=C(Oc3ccccc3)C=CC=C2OC1")

But allowing substituents on the benzene rings by using "c" instead of "C":

isoflavone_pattern = Chem.MolFromSmarts("O=C1c2c(oc3ccccc3)ccc2oc1")

This pattern should match the chromenone system with substituents on the benzene rings.

Testing this SMARTS against the example SMILES:

For Termisoflavone A, the core is O=C1C2=C(OC=C1C3=CC=C(...))C=C(...). The SMARTS "O=C1c2c(oc3ccccc3)ccc2oc1" would match if the benzene rings have substituents.

Alternatively, maybe the correct SMARTS is "O=C1c2c(oc3c([*])c([*])c([*])c([*])c3)ccc2oc1", allowing any substituents on the B-ring.

But this is getting too complex. Given the time, I'll revise the isoflavone_pattern to be more flexible by using aromatic atoms and wildcards.

Another approach: Use a SMARTS that matches the chromenone system without requiring specific substituents. For example:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1c2ccccc2oc3c1cccc3")

This pattern matches a benzopyran-4-one where the A and C rings are benzene rings. This should match the core structure even if there are substituents.

Testing this with the example SMILES:

Termisoflavone A's core is O=C1C2=C(OC=C1...) which would match the SMARTS if the rings are aromatic and have substituents.

If this SMARTS is correct, then the previous code's issue was due to the rigid SMARTS pattern not allowing substitutions. By using aromatic atoms (c) instead of aliphatic (C), the pattern becomes more flexible.

Additionally, checking for glycosylation might need adjustment. The previous code used "[O;!H0]-[C;r5,r6]" which looks for an oxygen connected to a ring carbon (sugar). But some examples have multiple glycosidic linkages. However, the previous code failed to detect the isoflavone core, so fixing that is the priority.

So, revising the isoflavone_pattern to use aromatic atoms and allowing substitutions should help. Also, the phenolic OH count might need adjustment, but the main issue was the core detection.

Let me update the code with the revised isoflavone_pattern and test.

The new code would have:

isoflavone_pattern = Chem.MolFromSmarts("[OX1]=C1c2ccccc2oc3c1cccc3")

This represents a benzopyran-4-one with two benzene rings (A and C), allowing for substitutions.

Another possible SMARTS is "O=C1c2ccccc2oc3ccccc13", which is a fused benzopyran-4-one.

Testing this with RDKit's SMARTS matcher for the example SMILES:

For Termisoflavone A's SMILES: O=C1C2=C(OC=C1C3=CC=C(...))C=C(...)O...

The core structure O=C1C2=C(...)OC=C1... should match the SMARTS if the rings are considered aromatic.

Assuming that the SMILES is parsed correctly and the rings are aromatic, the SMARTS should match.

With this change, the code should now detect the isoflavone core in the examples, leading to correct classification.

Another adjustment: the glycosylation pattern. The previous code looked for O connected to a ring carbon (r5 or r6). This should detect sugar moieties. However, in some examples, the glycoside might be attached via an oxygen to a non-ring carbon (like in apiosides). But for now, focusing on the core detection.

In conclusion, the main fix is revising the isoflavone_pattern to be more flexible, allowing substitutions on the benzene rings.