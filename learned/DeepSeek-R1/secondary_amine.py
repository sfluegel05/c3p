"""
Classifies: CHEBI:32863 secondary amine
"""
</think>

Okay, let's see why the previous code didn't work well. The main issue was a lot of false positives and false negatives. Let me think about why that happened.

Looking at the false positives, like LSM-12047 and others, the code incorrectly identified them as secondary amines. The code checks for nitrogen atoms with exactly three bonds and two carbon neighbors. But wait, secondary amines can have more than three bonds if the nitrogen is part of a conjugated system or has a double bond. For example, in an aromatic ring like pyridine, the nitrogen has two bonds but isn't a secondary amine. Also, the code didn't check if the nitrogen is in an amine group specifically. Some nitrogens might be in amides, nitro groups, or other functional groups which aren't amines.

The false negatives, like dobutamine, were missed because maybe the code didn't account for different bonding scenarios. For example, if the nitrogen has two carbons and one hydrogen but is part of a ring or has a double bond, the previous code might have missed it. Also, the code checks for degree 3, but in some cases, the nitrogen might have a different valence due to resonance.

So, the key points to fix are:
1. Ensure the nitrogen is actually in an amine group (not part of an amide, nitro, etc.).
2. Check the number of carbon neighbors (exactly two) and hydrogen count (at least one, but maybe not necessary if other substituents are present).
3. Consider the valence and possible double bonds that might affect the bonding structure.
4. Exclude nitrogens that are part of aromatic systems or other non-amine groups.

How to implement this in RDKit? Maybe use SMARTS patterns to specifically match the secondary amine structure. A secondary amine has a nitrogen connected to two carbons and one other group (could be another carbon or hydrogen). But the definition says two hydrocarbyl groups, which are carbon-based. So the nitrogen should have two carbon substituents and one hydrogen (or another substituent, but according to the definition, it's derived from ammonia by replacing two hydrogens with hydrocarbyl groups, so the third substituent would be a hydrogen).

Wait, but in reality, secondary amines can have three substituents (two carbons and one hydrogen) or maybe more if there's a charge (like in ammonium salts). But for the neutral case, it's two carbons and one hydrogen. However, in some structures, the hydrogen might not be present (e.g., if the nitrogen is part of a ring with two carbons and another substituent). Hmm, maybe the definition allows for any two hydrocarbyl groups and one other group (could be another hydrocarbyl, but that would make it tertiary). Wait, the definition says replacing two hydrogens with hydrocarbyl groups. So the remaining substituent on the nitrogen should be a hydrogen. But in practice, some secondary amines might have substituents other than hydrogen (like in rings), but according to the strict definition, maybe not.

Alternatively, maybe the correct approach is to look for a nitrogen with exactly two single bonds to carbon atoms and one hydrogen. But how to check for the hydrogen? Because in SMILES, explicit hydrogens might not be present. So using a SMARTS pattern that matches [NH](C)(C), where the nitrogen has two carbons and one hydrogen. But RDKit's SMARTS would need to account for implicit hydrogens.

So the SMARTS pattern for a secondary amine could be [N;H1;D3](C)(C), which means a nitrogen with three bonds (degree 3), one hydrogen, and two carbon neighbors. But maybe that's too restrictive because some secondary amines might have substituents other than hydrogen. Wait, according to the definition, the third substituent should be a hydrogen (since two hydrogens are replaced by hydrocarbyl). So the correct SMARTS would be [N;H1;X3]([#6])([#6]).

But wait, in some cases, the hydrogen might not be explicit. So using the SMARTS pattern [N;D3]([#6])([#6]) and then ensuring that the valence allows for a hydrogen. Alternatively, using the SMARTS pattern [N;H1;X3]([#6])([#6]).

Testing this pattern on the examples: For N-methylaniline (CNc1ccccc1), the nitrogen is connected to two carbons (methyl and benzene) and one hydrogen. So the SMARTS should match. For methamphetamine (CN[C@@H](C)Cc1ccccc1), the nitrogen is connected to two carbons (methyl and the chiral center) and a hydrogen. That should match.

But in the previous code, the check was for degree 3 and two carbon neighbors. But this would also match nitrogens in amides (like CON(C)(...)), but amides have a double bond to oxygen. So the code needs to exclude nitrogens that are part of other functional groups.

So the plan is:
1. Use a SMARTS pattern that matches the secondary amine structure: a nitrogen with exactly two carbon neighbors and one hydrogen (or the valence equivalent).
2. Exclude nitrogens that are part of other groups like amides, nitro, etc.

How to implement this in RDKit? Let's create a SMARTS pattern that matches [N]([#6])([#6]) where the nitrogen is not part of an amide, aromatic ring, etc.

The correct SMARTS pattern might be [N;!$(N-[O,S]);!$(N=*);!$(N#*);!$(N-C=[O,S,N])]([#6])([#6])[#0]. But I'm not sure. Alternatively, use [N;H1;X3]([#6])([#6]).

Wait, let's test the SMARTS pattern for secondary amines. The correct pattern should match NH(C)(C) but not N(C)(C)(C) (tertiary), or N in amides (like NC(=O)...).

So the SMARTS pattern [N;H1;X3]([#6])([#6]) should work. Let's see:

For N-methylaniline (CNc1ccccc1), the nitrogen has two carbons (methyl and benzene) and one H. So the pattern matches.

For an amide like CON(C)(...), the nitrogen is connected to two carbons and an oxygen (from the carbonyl), so the SMARTS pattern [N;H1;X3]([#6])([#6]) would not match because the degree is 3 but one neighbor is oxygen. So that's good.

But in the case of N-nitrosopiperidine (O=NN1CCCCC1), the nitrogen in the nitroso group (N-O) would have a double bond. The SMARTS pattern [N;H1;X3] requires three single bonds. So that would not match, but in the previous code, it was considered a false positive. Wait, in the error message, N-nitrosopiperidine was a true positive. But according to the definition, it's a secondary amine. Wait, the SMILES is O=NN1CCCCC1. The nitroso group is O=N-N. The nitrogen adjacent to the oxygen has a double bond. The other nitrogen (connected to the cyclohexane) has three bonds: one to the nitroso N, and two to carbons from the ring. Wait, but that nitrogen would have three bonds: two to carbons (from the ring) and one to the nitroso N. So in that case, the nitrogen in the ring has two carbons and one N neighbor. But according to the definition, a secondary amine should have two hydrocarbyl groups. The nitroso group's N is not a hydrocarbyl. So N-nitrosopiperidine should NOT be a secondary amine. But in the previous test, it was considered a true positive. Hmm, that's conflicting. So maybe the initial code had a mistake in considering that as a secondary amine, leading to a false positive.

So the problem with the previous code was that it considered any nitrogen with two carbons and three bonds as a secondary amine, regardless of other substituents. For example, in N-nitrosopiperidine, the nitrogen has two carbons (from the ring) and one N (from the nitroso group), so the code counted it as a secondary amine, but according to the definition, it's not because the substituents are not both hydrocarbyl.

Therefore, the code needs to ensure that both substituents on the nitrogen are carbons (hydrocarbyl groups), and the third substituent is a hydrogen (or possibly another group, but according to the definition, only two hydrogens are replaced by hydrocarbyl, leaving one hydrogen).

So the correct approach is to look for a nitrogen atom that is bonded to exactly two carbon atoms and one hydrogen. But since hydrogens are often implicit in SMILES, the code needs to account for that.

Using RDKit's SMARTS, the pattern [N;H1;D3]([#6])([#6]) would match a nitrogen with three bonds (degree 3), one implicit hydrogen, and two carbon neighbors. This should correctly identify secondary amines.

Additionally, we need to exclude nitrogens that are part of other functional groups, like amides (NC=O), enamines (N connected to double bond), aromatic nitrogens (like in pyrrole), etc.

So the revised code should:

1. Use a SMARTS pattern that matches the secondary amine structure: [N;H1;D3]([#6])([#6]).
2. Check if any such nitrogen exists in the molecule.
3. Exclude nitrogens that are part of other functional groups.

But how to exclude other functional groups? For example, amides have N connected to a carbonyl group. So the SMARTS pattern should exclude nitrogens adjacent to electronegative atoms or double bonds.

Alternatively, use a more specific SMARTS pattern that ensures the nitrogen is only connected to carbons and hydrogen. Let's refine the SMARTS:

The correct SMARTS pattern would be [N;H1;X3]([#6])([#6]), which matches a nitrogen with three bonds (valence 3), one hydrogen, and two carbon neighbors. This should exclude nitrogens in amides, nitro groups, etc.

Testing this pattern on the examples:

- N-methylaniline (CNc1ccccc1): The nitrogen has two carbons (methyl and benzene) and one H. Matches.
- Methamphetamine (CN[C@@H](C)Cc1ccccc1): The nitrogen has two carbons (methyl and the chiral C) and one H. Matches.
- N-nitrosopiperidine (O=NN1CCCCC1): The nitrogen in the ring has two carbons and one N (from nitroso). Does not match the SMARTS because the third neighbor is N, not H. So the code would correctly exclude it.

So using this SMARTS pattern should reduce false positives.

Now, implement this in code using RDKit's HasSubstructMatch method with the SMARTS pattern.

Another thing to consider: the SMARTS pattern should match any nitrogen that has exactly two carbon neighbors and one hydrogen, regardless of other substituents. Wait, but the definition requires exactly two hydrocarbyl groups (which are carbons), so the third substituent must be a hydrogen. So the SMARTS pattern [N;H1;D3]([#6])([#6]) should work.

Let me adjust the code to use this SMARTS pattern. Also, check for any matches in the molecule. If at least one such nitrogen exists, return True.

Testing this approach on the false positives and negatives:

False positives like N-nitrosopiperidine would not match the SMARTS pattern, so they would be excluded. False negatives like dobutamine would be matched if their nitrogen has two carbons and a hydrogen.

So the revised code would be:

- Parse the SMILES.
- Check for the presence of the SMARTS pattern.
- Return True if found, else False.

This should improve both false positives and negatives.

Also, need to handle cases where the hydrogen is explicit. But RDKit's SMARTS handles implicit hydrogens, so [N;H1] would match nitrogens with one implicit H.

Therefore, the code should use this SMARTS pattern and check for matches. Let's implement that.