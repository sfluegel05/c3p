"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
</think>

Okay, let's see why the previous code didn't work well. The main issue is that the SMARTS patterns weren't correctly identifying the beta-D-glucoside structure. The false negatives indicate that the patterns didn't match the actual examples provided. Let me think about what defines a beta-D-glucoside.

First, the glucose part: D-glucose has a specific stereochemistry. The hydroxyl groups on C2, C3, and C4 are in the same orientation. The anomeric carbon (C1) in beta configuration has the hydroxyl (or glycosidic oxygen) in the axial position. So the SMARTS pattern needs to capture the correct stereochemistry around the ring.

Looking at the previous SMARTS pattern: "[OX2;!H0]-[C@@H]1[C@@H]([OX2])[C@H]([OX2])[C@@H]([OX2])[C@H](C[OX2])O1". Wait, maybe the stereodescriptors are incorrect. Let me draw the beta-D-glucopyranose structure. The hydroxyls on C1 (beta) is up (axial), C2 is down (equatorial), C3 is up, C4 is down. So the chiral centers should have specific configurations.

Maybe the SMARTS pattern's stereodescriptors are not accurate. Let me adjust them. For beta-D-glucoside, the anomeric oxygen is connected to C1 with a specific configuration. The correct order of chiral centers around the ring (C1 to C5) in D-glucose is R, S, R, S. Wait, maybe using the @@ and @ symbols correctly is key here.

Alternative approach: Use RDKit's built-in function to find the anomeric carbon and check its configuration. But how? Maybe find the glycosidic oxygen and check the stereochemistry of the adjacent carbon (C1).

Another idea: Break down the problem into steps. First, identify all glucose residues in the molecule. Then check if any of them have the beta configuration at the anomeric position.

But how to detect glucose residues? Maybe look for a six-membered ring with five hydroxyl groups (or substituted) in the D-glucose configuration. Then check the anomeric position's stereochemistry.

Alternatively, use the SMILES notation for beta-D-glucopyranose: the anomeric oxygen is in the beta position. In SMILES, for beta, the hydroxyl (or OR group) on C1 is in the same direction as the C5 hydroxyl when drawn in the chair conformation. So in SMARTS, the C1 should have a specific configuration.

Wait, perhaps the correct SMARTS for beta-D-glucoside is "[C@H]1(O[C@@H]([C@@H]([C@H]([C@@H]1O)O)O)CO)". Let me check:

- The anomeric carbon (C1) is connected to the oxygen (glycosidic bond) with @ configuration.
- The C2 is @@ (since in D-glucose, the hydroxyl is on the opposite side from C1 in beta configuration).
Hmm, maybe I need to construct the SMARTS pattern carefully.

Another approach: Use the InChI code to identify beta-D-glucoside. But RDKit doesn't generate InChI directly. Alternatively, generate the canonical SMILES for beta-D-glucopyranose and compare substructures. But that might not work with substituents.

Wait, let's look at the example SMILES. For instance, genistein 7-O-beta-D-glucoside: OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)occ(-c2ccc(O)cc2)c3=O)[C@H](O)[C@@H](O)[C@@H]1O. The glucose part is O[C@@H](...)[C@H](O)... which corresponds to the beta configuration.

The key is that the glycosidic oxygen (O connecting to the aglycone) is on C1, and the stereochemistry around C1 and C2 is such that it's beta. The correct SMARTS should capture the C1 configuration.

Maybe the correct SMARTS for the beta-D-glucoside core is "[O;!H0]-[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1". Let me test this. The C1 (attached to the glycosidic O) has @@H (since in beta configuration, the O is axial, so the configuration would be R? Or S? Maybe I need to draw it out.

Alternatively, use a SMARTS that ensures the correct relative stereochemistry between C1 and C2. In beta-D-glucoside, the C1 and C2 hydroxyls are trans. So in the ring, the configuration around C1 and C2 should be such that their substituents are on opposite sides.

Maybe the correct SMARTS is "[OX2]-[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1". Let's see: the glycosidic O is connected to C1 (which is the anomeric carbon). The C1 has a certain configuration, and the adjacent C2 has another.

Alternatively, use the SMILES for beta-D-glucopyranose: O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1. The key is the @@ and @ in the positions.

So the SMARTS pattern should be "[OX2;!H0]-[C@@H]1[C@@H]([OX2;H0])[C@H]([OX2;H0])[C@@H]([OX2;H0])[C@H](CO)O1". This pattern looks for the anomeric oxygen (not a hydroxyl, since it's part of the glycosidic bond), then the correct stereochemistry around the ring.

But in the examples, some glucosides have substituents on the hydroxyls (like acetylations), so the hydroxyls might be replaced. So the pattern should allow any substituent on the oxygen atoms except the anomeric one. Wait, no: the anomeric oxygen is part of the glycosidic bond (O-R), while the other hydroxyls might be substituted. So the SMARTS should check for the anomeric oxygen (OX2, connected to two atoms, like in OR group), and the other positions can have any substituents (like OAc, etc.).

So the SMARTS pattern should be:

- An oxygen (OX2) connected to C1 (the anomeric carbon), which is part of a six-membered ring.
- The C1 has a specific stereochemistry (beta), which in SMARTS terms would be [C@@H] for the ring positions.
- The other carbons (C2, C3, C4, C5) have the correct D-glucose stereochemistry.

Putting this together, the correct SMARTS pattern should capture the ring structure with the correct stereochemistry. For beta-D-glucoside, the anomeric oxygen is in the axial position, which in the chair conformation of D-glucose would correspond to the C1 hydroxyl being trans to the C5 hydroxymethyl group.

But translating this into SMARTS requires precise stereochemical notation. Maybe using the SMARTS for beta-D-glucopyranose with the glycosidic oxygen:

The core structure is O-C1-C2-C3-C4-C5-O (ring), with specific stereochemistry.

Alternatively, use the following SMARTS pattern that captures the beta configuration:

"[OX2;!H0]-[C@@H]1[C@@H]([!#1])[C@H]([!#1])[C@@H]([!#1])[C@H](CO)O1"

Here, the anomeric oxygen (OX2;!H0) is connected to C1 with @@H (beta configuration). The other positions (C2, C3, C4) can have any substituents except hydrogen (since they are part of the glucose structure, which may have substituents). C5 is connected to the hydroxymethyl group (CO).

This pattern allows substituents on the hydroxyl groups (like acetylations) while ensuring the core glucose stereochemistry is correct.

Testing this pattern against the example SMILES:

For genistein 7-O-beta-D-glucoside (SMILES: OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)occ(-c2ccc(O)cc2)c3=O)[C@H](O)[C@@H](O)[C@@H]1O), the glucose part is O[C@@H](...), which should match the [C@@H] in the SMARTS. The other positions (C2, C3, C4) have hydroxyls or substituents, which the pattern allows.

Another example: cycasin (SMILES: C\[N+]([O-])=N/CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O). The glucose part is O[C@@H]1... which should match the pattern.

So, adjusting the SMARTS pattern to this more accurate version should capture the beta-D-glucoside structure. Also, ensuring that the D-configuration is checked by verifying the stereochemistry of the adjacent carbons (C2, C3, C4, C5).

Another improvement: Check that the glucose is in the D-form. In D-glucose, the C5 hydroxymethyl group is in a specific orientation. So the C5 should have the correct configuration. The SMARTS should include the correct stereochemistry around C5.

Wait, in D-glucose, the C5 (the last chiral center before the hydroxymethyl) has R configuration. So in the SMARTS, the C5 (which is part of the ring in pyranose form) should have a specific stereodescriptor.

Putting all together, the SMARTS pattern needs to accurately represent the stereochemistry of all chiral centers in the glucose ring. Let me construct it step by step.

The beta-D-glucopyranose structure in SMILES is O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1. The anomeric carbon (C1) has @@H, C2 has @H (since in D-glucose, the hydroxyls are on alternating positions). Wait, let me confirm:

In the beta-D-glucopyranose ring (chair conformation), the hydroxyls are arranged as follows: C1 (axial, beta), C2 (equatorial), C3 (axial), C4 (equatorial), and C5 (axial, leading to the hydroxymethyl group). The chiral centers (C2-C5) have specific configurations.

The correct SMARTS pattern should enforce the D-glucose stereochemistry. So the pattern would be:

"[OX2;!H0]-[C@@H]1[C@H]([!#1])[C@@H]([!#1])[C@H]([!#1])[C@@H](CO)O1"

Breaking it down:

- OX2;!H0: The glycosidic oxygen (not a hydroxyl).
- [C@@H]1: C1 with @@ configuration (beta).
- [C@H]([!#1]): C2 with @ configuration.
- [C@@H]([!#1]): C3 with @@ configuration.
- [C@H]([!#1]): C4 with @ configuration.
- [C@@H](CO): C5 with @@ configuration, connected to CH2OH (D-configuration).

This pattern ensures the correct stereochemistry for D-glucose in the beta configuration.

Testing this against the example SMILES should now correctly identify beta-D-glucosides. Additionally, this pattern allows substituents on the hydroxyl groups (since [!#1] matches any non-hydrogen atom, like OAc groups).

Another consideration: The hydroxymethyl group (C5) should be present. The pattern includes [C@@H](CO), ensuring that C5 is connected to CH2OH (or substituted, as in some examples).

By refining the SMARTS pattern to accurately capture the stereochemistry of all chiral centers in the beta-D-glucoside, the previous false negatives (like genistein glucoside, cycasin) should be correctly identified. The false positives from the previous attempt, which were other glycosides with different configurations (like alpha-linked or non-glucose sugars), would now be excluded because the pattern enforces the specific D-glucose stereochemistry.

Finally, to handle large molecules with multiple rings, the function should iterate through all matches to ensure at least one matches the beta-D-glucoside pattern. Also, check that the glucose is indeed present (like having six oxygen atoms in the ring, but substituents may vary).

In conclusion, the key improvement is using a more precise SMARTS pattern that accurately represents the stereochemistry of beta-D-glucoside, including all chiral centers in the glucose ring. This should reduce false negatives and prevent false positives from other glycoside types.