"""
Classifies: CHEBI:33709 amino acid
"""
This pattern should match any amine (primary, secondary, tertiary) that is not part of an amide.

So the steps are:

1. Check that the molecule has at least one carboxylic acid group.
2. Check that the molecule has at least one amine group (not part of an amide).

Now, implementing this in RDKit:

First, parse the SMILES. Then:

- Check for carboxylic acid: using the SMARTS "[CX3](=O)[OX2H1]". This matches COOH groups.
- Check for amine groups using the SMARTS pattern that excludes amides.

Additionally, the molecule should not be something like a beta-alanine (NH2-CH2-CH2-COOH), which is an amino acid. Beta-alanine's SMILES is NCC(=O)O. Wait, no: beta-alanine is NH2-CH2-CH2-COOH, which is NCCCC(=O)O. Wait, no: beta-alanine is NH2-CH2-CH2-COOH, SMILES: NCCCC(=O)O. So in this case, the amine is a primary amine (NH2) connected to a CH2-CH2-COOH chain. So the amine is not part of an amide, so it should be detected.

Another example: 3-aminopropionic acid (beta-alanine) should be detected as an amino acid.

Now, let's test the code with some examples.

For glycine (SMILES: NCC(=O)O), the code should find a carboxylic acid and an amine group. So returns True.

For acetylglycine (SMILES: CC(=O)NCC(=O)O), there is a carboxylic acid (the terminal COOH), but the amine is part of an amide (the acetyl group). So the code should not find any amine groups (other than the amide), so returns False.

For N-acetyl-L-lysine (SMILES: CC(=O)NCCCC[C@H](N)C(=O)O), there's a carboxylic acid (C(=O)O) and an amine group (the [C@H](N) part), which is a primary amine not part of an amide. So returns True.

Another example: 3-aminobenzoic acid (SMILES: Nc1cccc(c1)C(O)=O). It has a carboxylic acid and an amino group on the benzene ring. So returns True.

Now, let's code this.

Potential issues:

- The carboxylic acid SMARTS might not match all cases. For example, if the carboxylic acid is deprotonated (like COO-), but in SMILES, it's written as COO. But the given examples use C(O)=O, which is correct. The SMARTS "[CX3](=O)[OH]" should match COOH groups.

Another possible edge case: molecules with multiple carboxylic acids and multiple amines. But the function only needs at least one of each.

So, the code steps:

1. Parse the SMILES into a molecule. If invalid, return False and reason.
2. Check for at least one carboxylic acid group using the SMARTS pattern.
3. Check for at least one amine group using the SMARTS pattern that excludes amides.
4. If both are present, return True with the reason. Otherwise, return False with the appropriate reason.

Now, implementing this in Python:

The code should import RDKit modules, define the function, and perform these checks.

Testing the previous code's error: the user's last attempt had an error due to incorrect indentation. So the code must be properly indented.

Putting it all together:

The function is_amino_acid takes a SMILES string, converts it to a molecule. If that fails, return False. Then check for carboxylic acid and amine groups.

So the code would be: