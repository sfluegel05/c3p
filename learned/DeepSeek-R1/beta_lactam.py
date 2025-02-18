"""
Classifies: CHEBI:35627 beta-lactam
"""
[ND1]=C-O might not be right. Wait, the amide is N connected to C=O. So in the ring, there's a bond N-C(=O). So the SMARTS could be "[nH0]C(=O)@[r4]" to find a nitrogen (with one H, like in a secondary amide) connected to a carbonyl carbon in a four-membered ring.

Alternatively, a SMARTS pattern that looks for a four-membered ring where there's an N and a C=O connected in the ring. For example:

"[nH0]1-C(=O)-*-1" where the ring is four members. But SMARTS notation requires specifying the exact connectivity.

Alternatively, use a SMARTS pattern that finds any four-membered ring where there's an N and a C=O in the ring, adjacent to each other.

After some research, the correct SMARTS pattern for a beta-lactam (four-membered lactam) is "[nH0]1C(=O)CC1". But wait, that's what the original code used. But in the case of fused rings, like penicillins, the nitrogen is part of a bicyclic system. Let's look at penicillin's structure: the lactam ring is part of a bicyclo[3.2.0] system. The nitrogen is in a four-membered ring (the lactam) and also part of a five-membered ring. So in the SMILES for penicillin, the nitrogen is connected to three atoms (two in the four-membered ring and one in the five-membered ring). The original SMARTS pattern [NH]1C(=O)CC1 would not match because the nitrogen has three bonds (so it's N with no H, i.e., [NH0] instead of [NH]).

Ah! The original code used "[NH]1C(=O)CC1", which requires the nitrogen to have a hydrogen. But in penicillin, the nitrogen is part of a bicyclic system and doesn't have a hydrogen. So the SMARTS pattern should use [NH0] instead of [NH]. That's a critical mistake. The correct SMARTS should account for the nitrogen's valence in the ring. So changing the SMARTS pattern to "[nH0]1-C(=O)-C-C-1" would match the nitrogen in the ring without requiring a hydrogen.

So the previous code's SMARTS was looking for [NH], which would miss cases where the nitrogen is part of a fused ring and has no hydrogen. Changing the nitrogen to [nH0] (which allows any nitrogen in a ring with no hydrogen) would fix that.

Let me test this. For example, in penicillin's SMILES: [H][C@@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(c1ccccc1)S(O)(=O)=O)C(O)=O. The lactam ring is part of the bicyclo system. The nitrogen in the lactam (N1) is connected to C(=O) and two other carbons (part of the rings). So the SMARTS pattern [nH0]1C(=O)CC1 would match because the nitrogen is in a four-membered ring with C(=O), C, C.

So the fix is to adjust the SMARTS pattern to allow the nitrogen to have no hydrogen and be part of a four-membered ring with the amide bond.

Additionally, some beta-lactams might have substituents on the ring atoms. The original SMARTS might not account for substituents. For example, the C adjacent to the carbonyl might have branches. So the SMARTS pattern should allow any atoms in the ring as long as the N and C=O are present.

Another possible issue is that the original code's SMARTS requires the exact sequence N-C(=O)-C-C, but in reality, the ring could have different connectivity. For example, the N could be connected to C-C(=O)-C. So a more flexible SMARTS pattern that allows any arrangement of N and C=O in a four-membered ring would be better.

To capture all possible four-membered rings containing an amide bond, regardless of the order, we can use a SMARTS pattern that checks for the presence of both a nitrogen and a carbonyl group within the same four-membered ring.

An improved approach would be:

1. Find all atoms in four-membered rings.
2. For each four-membered ring, check if it contains a nitrogen and a carbonyl carbon (C=O) that are adjacent.

But implementing this with SMARTS might be complex. Alternatively, use a SMARTS pattern that looks for any four-membered ring where a nitrogen is adjacent to a carbonyl carbon.

The SMARTS pattern could be "[nH0]C(=O)@[r4]" which means a nitrogen (without H) connected to a carbonyl carbon, both in a four-membered ring. However, the exact syntax might need adjustment.

Alternatively, use a SMARTS pattern that looks for a four-membered ring with a nitrogen and a carbonyl group connected in any order. For example:

"[r4]@[nH0]@C(=O)" — but SMARTS requires specifying the exact connectivity.

Alternatively, use a SMARTS pattern that allows any four-membered ring containing both an N (without H) and a C=O, connected in the ring.

After some research, the correct SMARTS pattern would be "[nH0;r4]C(=O);@[r4]". But I'm not sure. Let me test this.

Another approach is to use the SMARTS "[nH0]1-C(=O)-*-1" where the * represents any two atoms to complete the four-membered ring. For example:

beta_lactam_pattern = Chem.MolFromSmarts("[nH0]1-C(=O)-C-C-1")

But this still requires the N and C=O to be adjacent. However, in some cases, the N and C=O might not be adjacent but still part of the same four-membered ring. Wait, no—the amide bond is N-C=O, so they must be adjacent. Because in a lactam, the amide is part of the ring, so the N is directly bonded to the carbonyl carbon. Therefore, the correct SMARTS should have N adjacent to C=O in the ring.

So the correct SMARTS pattern should be "[nH0]C(=O)@[r4]" — but how to write that.

Alternatively, the SMARTS pattern for a four-membered ring where the N and C=O are adjacent:

"[nH0]1-C(=O)-*-*-1"

But in RDKit's SMARTS, this would be written as "[nH0]1-C(=O)-*-*-1", but the * might not be necessary. Let's try:

beta_lactam_pattern = Chem.MolFromSmarts("[nH0]1-C(=O)-*-*-1")

But I'm not sure. Alternatively, using explicit atoms:

"[nH0]1-C(=O)-C-C-1"

But this requires the N to be connected to C(=O), then two more carbons. However, in some beta-lactams, there might be other atoms (like sulfur in penicillins). Wait, no—the lactam ring itself is four-membered with N, C=O, and two other carbons. Even in penicillins, the lactam ring is part of the bicyclo system, but the lactam ring itself is four atoms: N, C=O, and two other carbons (possibly with substituents).

Wait, in penicillin, the lactam ring is part of a bicyclo[3.2.0] system. The four-membered ring would be N connected to C=O, then two more carbons connected back to the N. So the SMARTS pattern "[nH0]1-C(=O)-C-C-1" should match that.

But in the SMILES for penicillin, the N is part of the bicyclo system. Let's take ampicillin's SMILES: [H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O. The N1 is part of the four-membered ring (positions 1, 2, 3, 4?), with C(=O) at position 2. So the SMARTS pattern should match.

So the problem in the previous code was using [NH] (which requires a hydrogen on the nitrogen) instead of [nH0] (which allows the nitrogen to have no hydrogen, as in fused rings).

So the corrected SMARTS should be "[nH0]1-C(=O)-C-C-1".

Let me adjust the code to use this SMARTS pattern.

Another thing to consider: the amide nitrogen might be part of a bridge in a bicyclic system. The SMARTS should still match as long as the four-membered ring contains the N and C=O.

Testing with the example that was a false negative: doripenem's SMILES has a part where the N is in a four-membered ring with C=O. Let's see:

Doripenem SMILES: C[C@@H](O)[C@@H]1[C@H]2[C@@H](C)C(S[C@@H]3CN[C@@H](C3)CNS(N)(=O)=O)=C(N2C1=O)C(O)=O

Looking at the structure, the lactam ring is part of a bicyclic system. The N is in a four-membered ring connected to C=O. So using the corrected SMARTS should match.

Now, the false positives: for example, Streptoverticillinone's SMILES is O=C1N[C@@](C(=O)O)(CC2=CC=C(O)C=C2)C1. This forms a four-membered ring with N and C=O. According to the definition, this should be a beta-lactam. But the test case marked it as a false positive. Wait, that's confusing. If the definition is any four-membered lactam, then this should be a beta-lactam. But perhaps the test cases are incorrect, or there's a mistake in the classification. Alternatively, maybe the test cases are considering only certain subclasses of beta-lactams (like antibiotics), but the code should follow the strict chemical definition.

Assuming the code should follow the strict definition, then the previous code's issue was the SMARTS pattern requiring [NH], which missed many cases. Fixing that would reduce false negatives but might include some true positives that were previously false positives. However, according to the given outcomes, some structures like Azetirelin and Streptoverticillinone were marked as false positives, but according to the definition, they are beta-lactams. This suggests that the test cases might have errors, but the code should adhere to the chemical definition.

Therefore, the main correction needed is to adjust the SMARTS pattern to [nH0]1C(=O)CC1 instead of [NH]1C(=O)CC1.

So the revised code would use the corrected SMARTS pattern and check for any four-membered ring containing an amide bond, regardless of the nitrogen's hydrogenation state.

Additionally, to handle cases where the lactam is part of a larger ring system, the code should use the SMARTS pattern that matches the four-membered ring with the amide, regardless of other substituents.

Let's implement this change and see if it improves the results.