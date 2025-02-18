"""
Classifies: CHEBI:33704 alpha-amino acid
"""
</think>

Okay, let's figure out why the previous code didn't work well. The main issues were false positives and false negatives. Let's look at the errors.

False positives happened with molecules that have amino and carboxylic groups but aren't alpha-amino acids. For example, dipeptides like Val-Phe were incorrectly classified. The original code checks if any carboxylic acid has an adjacent carbon with an amino group, but in dipeptides, the peptide bond connects two amino acids. The code might be detecting the amino group from the second amino acid attached to the carboxylic group of the first, leading to a false positive.

False negatives occurred when the amino group wasn't detected properly. Cases like fusarinine zwitterion might have the amino group in a different form (like a zwitterion) that the original pattern didn't catch. Also, deuterated compounds like glycine-d5 might have hydrogen isotopes that affected the hydrogen count check.

So, the main problems are:
1. The code checks for any carboxylic acid with an adjacent amino group, which can include peptides.
2. The amino group detection might miss some variations (like zwitterions or substituted amines).
3. Not ensuring that the amino and carboxylic groups are on the same carbon chain (like in branched structures).

To fix this:
- Need to ensure that the amino group is directly attached to the alpha carbon, which is adjacent to the carboxylic acid's carbon.
- Exclude cases where the amino group is part of an amide (like in peptides) by checking for adjacent carbonyl groups.
- Consider different protonation states (like zwitterions) by looking for any nitrogen connected to the alpha carbon, regardless of hydrogen count.
- Verify that the carboxylic acid and amino group are part of the same backbone, not in separate branches.

Let's adjust the approach:
1. Use SMARTS patterns to find the alpha carbon (connected to both the carboxylic acid and the amino group).
2. Check that the amino group isn't part of an amide or other conjugated system.
3. Ensure the carboxylic acid is present and properly detected, even in deprotonated forms.

Revised steps:
- Find all carboxylic acid groups (including deprotonated forms).
- For each carboxylic acid, check the adjacent carbon (alpha carbon) for an amino group.
- The amino group should be a primary, secondary, or tertiary amine but not part of an amide.
- Exclude structures where the amino group is part of a peptide bond (amide).

Now, let's implement this with RDKit SMARTS patterns. The key SMARTS would be something like [C](=O)[O] connected to a carbon that has a nitrogen attached, but that nitrogen shouldn't be part of an amide.

Maybe a SMARTS pattern like:
- Carboxylic acid: [CX3](=O)[OX2H0,O-]
- Alpha carbon connected to carboxylic acid and amino group: [CX4H][CX3](=O)[OX2H0,O-]

Wait, perhaps a better approach is to create a SMARTS that matches the alpha-amino acid structure directly. For example:

The core structure is:
- Carboxylic acid group: C(=O)O
- Adjacent carbon (alpha carbon) connected to an amino group (NH2, NHR, NR2)

So the SMARTS could be: [NH0;H2,H1;!$(NC=O)][CH2][C](=O)[OH0-,O]

Breaking it down:
- [CH2] is the alpha carbon (connected to the carboxylic carbon and the amino group).
- [C](=O)[OH0-,O] is the carboxylic acid (including deprotonated forms).
- [NH0;H2,H1;!$(NC=O)] is the amino group (not part of an amide). The !$(NC=O) ensures the nitrogen isn't adjacent to a carbonyl (amide).

This should match the alpha-amino acid structure while excluding amides.

Testing this SMARTS pattern on the examples should help. Let's adjust the code to use this pattern and see if it reduces false positives and negatives.

Another consideration: the original code loops through each carboxylic acid and checks neighbors. But using a single SMARTS pattern that captures the entire required structure might be more efficient and accurate.

So, the plan is:
1. Use a SMARTS pattern that matches the alpha-amino acid structure.
2. Check for at least one match of this pattern in the molecule.
3. Also, ensure there's a carboxylic acid group (in case the SMARTS misses some variants).
4. Exclude amide bonds by ensuring the nitrogen isn't part of an amide.

Let's implement this and test against the given examples.