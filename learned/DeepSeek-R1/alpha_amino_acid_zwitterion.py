"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
</think>

Okay, let's see why the previous code didn't work well. The main issue was false positives and negatives. Looking at the errors:

False positives happened when the code found [NH3+] and [O-]C(=O) on adjacent carbons but the structure wasn't an alpha-amino acid. For example, D-serine was correctly identified as a zwitterion but is an alpha-amino acid. Wait, but D-serine should be a valid zwitterion. Hmm, maybe some false positives were actually correct. But the problem mentions some entries were wrongly classified. Let me check the examples again.

Wait, in the false positives list, there's "D-serine zwitterion" marked as wrong. But according to the definition, alpha-amino acid zwitterions include both D and L forms. So maybe the code should allow stereochemistry but the test case considered it wrong. Or perhaps the definition requires the amino group to be on the alpha carbon relative to the carboxylate, which the code checks.

Looking at the false negatives: trans-4-hydroxy-L-proline zwitterion was missed because it has [NH2+] instead of [NH3+]. The original code only looked for [NH3+], but proline's amino group is part of a ring, leading to [NH2+]. So the code needs to account for secondary amines like in proline derivatives.

Another false negative was (2S,3R)-diazaniumylbutanoate, where the ammonium is not on the alpha carbon. Wait, but that's a different structure. The code requires the amino and carboxyl to be on adjacent carbons. So maybe that's correctly a false negative.

So the main issues are:

1. Missing proline-like zwitterions where the ammonium is [NH2+] (secondary amine).
2. Incorrectly allowing cases where the amino and carboxyl are adjacent but not in the alpha position (but how? The code checks adjacent carbons, which is the alpha position).
3. Maybe some stereochemistry issues, but the problem's test cases might have errors.

To fix the proline issue, the code should look for both [NH3+] and cyclic [NH2+] groups that are part of the alpha carbon. So adjust the ammonium pattern to include secondary amines in a ring.

Also, in the false positives, some molecules might have the ammonium and carboxylate groups adjacent but not in the correct position. For example, if the carboxylate is not the main one. Wait, the definition says "transfer of a proton from the carboxy to the amino group of any alpha-amino acid". So the molecule must be derived from an alpha-amino acid, meaning the amino and carboxyl groups are on the same carbon (the alpha carbon). Wait, no: in alpha-amino acids, the amino group is on the carbon adjacent to the carboxyl. So the code's original logic of checking adjacent carbons is correct.

Wait, alpha-amino acids have the amino group attached to the carbon adjacent to the carboxyl group (the alpha carbon). So the code checks that the ammonium's carbon (alpha) is adjacent to the carboxyl carbon. That's correct.

But in proline, the amino group is part of a ring. So the ammonium would be [NH2+] because it's a secondary amine. The original code misses this because it only looks for [NH3+].

So to handle proline zwitterions, the code should also check for [NH2+] in a ring. But how to differentiate between a secondary amine in a zwitterion vs other cases?

Another approach: check that the nitrogen with positive charge (whether [NH3+] or [NH2+]) is bonded to the alpha carbon, which is adjacent to the carboxyl carbon.

So steps to improve:

1. Update the ammonium pattern to include both [NH3+] and [NH2+] where the nitrogen is positively charged and part of the amino group attached to the alpha carbon.

Wait, but how to represent [NH2+] in a ring? The SMARTS pattern for a secondary amine with positive charge would be [NH2+].

But in proline, the nitrogen is part of a five-membered ring. So the code should find [NH2+] that's bonded to two carbons (part of the ring) and has a positive charge.

So modify the ammonium pattern to look for any positively charged amino group attached to a carbon.

Wait, the definition is about zwitterions from alpha-amino acids. In such cases, the amino group (now protonated) is on the alpha carbon. So the code needs to ensure that the nitrogen is bonded to the alpha carbon, which is adjacent to the carboxyl carbon.

So perhaps the problem with the previous code was that it only looked for [NH3+], missing cases where the nitrogen is [NH2+] (like proline).

So adjusting the ammonium pattern to include [NH2+] when it's part of the alpha carbon.

Let me adjust the code:

Change the ammonium_pattern to [N+;H2,H3].

In SMARTS, [N+;H2,H3] would match nitrogen with positive charge and 2 or 3 hydrogens (like NH3+ or NH2+ where one H is replaced by a bond in a ring).

Wait, but in proline's zwitterion, the nitrogen has two bonds to carbons (in the ring) and one H, leading to a positive charge. So the SMARTS would be [N+;H1], but that might not be general.

Alternatively, use [N+]([C])([C])[C] for secondary amines in a ring, but that's too specific.

Alternatively, look for any nitrogen with a positive charge that is bonded to a carbon (the alpha carbon), which in turn is adjacent to the carboxyl carbon.

So the steps would be:

1. Find any nitrogen with a positive charge (regardless of H count) that is bonded to a carbon (alpha carbon).
2. Check that this alpha carbon is adjacent to the carboxyl carbon.

This approach would capture both primary ([NH3+]) and secondary ([NH2+]) amines, as long as they are attached to the alpha carbon next to the carboxyl.

So modifying the code:

- Change the ammonium_pattern to [N+], which matches any positively charged nitrogen.
- Then, for each such nitrogen, check if it's bonded to a carbon (alpha carbon), which is adjacent to the carboxyl carbon.

This should include proline derivatives where the nitrogen is part of a ring.

Also, need to ensure that the carboxylate is the main one (i.e., the one adjacent to the alpha carbon), not another carboxylate in the molecule.

Another issue in the previous code: it returns True as soon as any pair of ammonium and carboxylate groups are adjacent. But the molecule might have multiple such groups. However, the definition requires that it's derived from an alpha-amino acid, which has exactly one amino and one carboxyl group. But some examples have multiple charges (like the dizwitterion), which the code might have handled.

Wait, the example includes "epsilon-(gamma-glutamyl)lysine dizwitterion" which has two ammonium and two carboxylate groups. So the code should check that at least one pair of ammonium and carboxylate are on adjacent carbons, but others can be elsewhere.

But according to the definition, the zwitterion is formed by transferring a proton from the carboxy to the amino group of any alpha-amino acid. So the core structure should have the alpha-amino acid part as zwitterion, even if there are other groups.

So the code's original approach of checking for at least one such pair is correct.

So, modifying the ammonium pattern to [N+] to capture all positively charged nitrogens, then checking if any of them are attached to a carbon (alpha) adjacent to a carboxyl carbon.

Testing this approach:

For proline zwitterion (trans-4-hydroxy-L-proline zwitterion), SMILES: O[C@H]1C[NH2+][C@@H](C1)C([O-])=O.

The nitrogen is [NH2+], part of a ring. The code would find this nitrogen ([N+]), check its attached carbon (the alpha carbon, which is part of the ring), and see if that carbon is adjacent to the carboxyl carbon.

In this case, the alpha carbon (attached to N) is connected to the carboxyl carbon (C=O). So the code should return True.

In the previous code, the ammonium_pattern was [NH3+], so it missed this case.

So changing the ammonium pattern to [N+] should fix this.

Another adjustment: in the previous code, the carboxylate pattern was [O-]C(=O), but sometimes the carboxylate could be written as C(=O)[O-]. The SMARTS pattern [O-]C(=O) should match either way, as SMARTS is substructure-based.

But let's confirm: the SMILES for carboxylate is C([O-])=O, which as a SMARTS would be [CX3](=[OX1])[O-]. Alternatively, using [O-]C(=O) would match O- attached to C=O.

Wait, the SMARTS [O-]C(=O) would match O- connected to a carbon that is double-bonded to O. So that should correctly identify the carboxylate group.

So the carboxylate pattern is okay.

Another possible issue: in some cases, the carboxylate might be part of an ester or other group, but according to the definition, it should be the main carboxyl group from the amino acid.

But the code can't distinguish that; it just checks for presence. However, the definition requires that the zwitterion is formed from the alpha-amino acid, so the carboxylate must be the original one adjacent to the amino group.

So the code's logic is correct in that aspect.

Now, let's adjust the code:

1. Change the ammonium pattern to [N+], which includes all positively charged nitrogens.
2. For each such nitrogen, find the attached alpha carbon (must be a carbon directly bonded to N).
3. Check if any of the alpha carbons are adjacent to a carboxyl carbon.

This should handle both primary and secondary amines.

Testing the false negatives:

- trans-4-hydroxy-L-proline zwitterion has [NH2+], which the new pattern would catch. The alpha carbon (attached to N) is part of the ring and adjacent to the carboxyl carbon. So the code should return True.

Another false negative was (2S,3R)-diazaniumylbutanoate, which has the ammonium and carboxylate on the same carbon? Let's look at its SMILES: O=C([O-])[C@@H](N)[C@@H](C)[NH3+]. Wait, the carboxylate is on one carbon, and the ammonium is on another carbon. The structure is NH3+ - C - C - COO-. So the alpha carbon (attached to NH3+) is adjacent to the next carbon, which is attached to the carboxylate. So the alpha carbon and carboxyl carbon are not adjacent. Hence, the code correctly returns False, which explains why it was a false negative but actually correct.

So the code's original logic was correct in that case, but the test case might have been wrong.

Another false negative was thioproline zwitterion (SMILES: C(C1[NH2+]CSC1)(=O)[O-]). Here, the ammonium is [NH2+] in a ring. The alpha carbon is the one attached to N, which is connected to the carboxyl carbon. So the code should detect this.

In the previous code, the ammonium pattern was [NH3+], so it missed this. With the new pattern [N+], it would find the [NH2+] and proceed.

So the modified code should fix these issues.

Another adjustment: in the previous code, when checking the alpha carbon, it looked for any carbon attached to the nitrogen. But in cases where the nitrogen is part of a ring (like proline), the alpha carbon is the one that's also adjacent to the carboxylate.

So the code's logic to find the alpha carbon as any carbon attached to N is correct.

Now, updating the code accordingly.

Testing the false positives:

For example, D-serine zwitterion ([NH3+][C@H](CO)C([O-])=O) was marked as a false positive, but according to the definition, it should be a valid alpha-amino acid zwitterion regardless of stereochemistry. So perhaps the test case was incorrect. The code correctly identifies it, but the test case wrongly considers it a false positive. However, the user's data says it's a false positive, which suggests that the code should not classify it. But according to the problem statement, the examples provided include L-lysine, D-tryptophan, etc. So stereochemistry shouldn't matter. Therefore, the test case might have an error, but the code should accept both D and L forms.

Alternatively, maybe the false positive was due to another reason. Looking at the false positives list, some entries have multiple charges or other groups. For example, gamma-Glu-Gly(1-) has an ammonium and carboxylate but maybe in a different arrangement.

The code's logic would check if any ammonium group is adjacent to a carboxylate. In gamma-Glu-Gly, the ammonium might be on the glutamate part, and the carboxylate on the glycine. Depending on the structure, they might not be adjacent. But the code would need to verify.

In any case, the main fix is to handle [NH2+] cases.

Another possible issue: the code returns True as soon as any pair is found, but some molecules might have multiple ammonium and carboxylate groups. However, the definition allows for this as long as at least one pair forms the alpha-amino acid zwitterion.

So the code's approach is acceptable.

In conclusion, the main fix is to adjust the ammonium pattern to include any positively charged nitrogen, not just [NH3+], and ensure that the alpha carbon is adjacent to the carboxyl carbon. This should capture proline derivatives and other secondary amine zwitterions.