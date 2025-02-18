"""
Classifies: CHEBI:46640 diketone
"""
</think>

Okay, let's figure out why the previous attempt failed. The error mentioned was "unexpected indent" which suggests there was a syntax issue in the code, like incorrect spacing or tabs. Looking at the example provided, the user pasted a long list of SMILES strings, which might have been part of the previous code instead of the actual function. The actual code for is_diketone wasn't shown, so maybe the previous answer included those SMILES strings by mistake, causing a syntax error.

Now, the task is to write a function called is_diketone that checks if a molecule has exactly two ketone groups. A diketone has two ketone functionalities, which are carbonyl groups (C=O) where the carbonyl carbon is connected to two other carbons (not part of a carboxylic acid or ester, etc.).

So, the plan is:

1. Parse the SMILES into a molecule.
2. Check for invalid SMILES.
3. Find all ketone groups. A ketone is a carbonyl group (C=O) where the carbon is connected to two other carbons. So, the SMARTS pattern would be [CX3](=O)[#6] where the carbon has exactly three bonds (double bond to O and two single bonds to other carbons).
4. Ensure there are exactly two such groups.
5. Make sure these ketones aren't part of other functional groups like acids or esters. For example, in a carboxylic acid (O=C-OH), the carbonyl is not a ketone. So, we need to exclude cases where the carbonyl is adjacent to an oxygen connected to a hydrogen (like in acids) or part of an ester.

Wait, but how to differentiate ketones from other carbonyl-containing groups? The SMARTS pattern [CX3](=O)[#6] would match ketones, but also amides or other groups where the carbonyl is connected to non-oxygen atoms. But maybe for simplicity, just count the number of ketone groups using that pattern and ensure there are exactly two.

Alternatively, use the SMARTS pattern that excludes adjacent oxygen. For example, ketone is [C;$(C(=O)[#6])](=O)[#6], but maybe a better approach is to use the RDKit's built-in functional group detection or ensure that the carbonyl is not part of an acid, ester, amide, etc.

Alternatively, the SMARTS for a ketone could be [CX3]=[OX1] where the carbon is connected to two other carbons. Wait, maybe the correct SMARTS is [CX3]([#6])([#6])=O. That would match a carbon connected to two carbons and a double-bonded oxygen.

Testing that: for acetone (CC(=O)C), the pattern [CX3]([#6])([#6])=O would match. For a carboxylic acid (C(=O)O), the carbonyl carbon is connected to an oxygen, so it wouldn't match.

So, the SMARTS pattern for a ketone group is [CX3H0](=[OX1])[#6]. The H0 ensures that the carbon has no hydrogens (since in ketones, the carbonyl carbon is typically R2C=O, so no H). But maybe some cyclic ketones might have H? Not sure. Alternatively, use [CX3](=[OX1])[#6] and ensure that the adjacent atoms are not oxygen (like in esters or acids).

Alternatively, find all carbonyl groups and then check if they are ketones by their environment.

But perhaps the simplest approach is to use the SMARTS pattern [CX3]=O where the carbon is connected to two other carbons. Wait, SMARTS is [CX3]([#6])([#6])=O. Let's test this:

For acetone (CC(=O)C), the central carbon is connected to two carbons and a double-bonded O. So it matches.

For propanal (CCC=O), the carbonyl carbon is connected to one carbon and two Hs. The SMARTS [CX3]([#6])([#6])=O would not match because the carbonyl carbon has only one carbon neighbor. So that's an aldehyde, not a ketone.

For a carboxylic acid (C(=O)O), the carbonyl carbon is connected to an O, so it wouldn't match.

So using the SMARTS pattern [CX3]([#6])([#6])=O would correctly identify ketones.

So, the steps are:

- Convert SMILES to molecule.
- Check validity.
- Find all substructures matching the ketone pattern.
- If exactly two matches, return True, else False.

But wait, some molecules might have more than two ketones. For example, triketones would have three. So the function should check that there are exactly two.

So the code would look something like:

from rdkit import Chem

def is_diketone(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    ketone_pattern = Chem.MolFromSmarts('[CX3]([#6])([#6])=O')
    matches = mol.GetSubstructMatches(ketone_pattern)
    if len(matches) == 2:
        return True, "Exactly two ketone groups found"
    else:
        return False, f"Found {len(matches)} ketone groups, expected 2"

But wait, maybe the pattern is too strict. For example, in a cyclic ketone like cyclohexanone, the carbonyl carbon is connected to two carbons (in the ring). So the pattern should match.

Testing the example SMILES given:

Take 2,5-hexanedione: CC(=O)CCC(C)=O. This has two ketone groups. The code should find two matches.

Another example: curcumin has two ketone groups (the two CC(=O) groups in the middle). The code should find two.

But wait, some of the examples have conjugated systems. Let's see.

Another example: acetylacetonate ([H][C-](C(C)=O)C(C)=O). Here, the two ketone groups are present. The pattern [CX3]([#6])([#6])=O would match each of the two carbonyl carbons. Each is connected to two carbons (the central C is connected to two CH3 groups and the O). So yes, two matches.

So this pattern should work.

But what about cases where the ketone is part of a ring? For example, bruceolline E (SMILES: CC1(C)C(=O)C(=O)c2c1[nH]c1ccccc21). Here, there are two ketone groups (the two =O). Each carbonyl carbon is connected to two carbons (from the ring and the adjacent groups). So the pattern would match both, so len(matches) is 2, which is correct.

So the code should work. However, in the previous attempt, the code might have had incorrect syntax, like the SMILES list being part of the code. So the correct approach is to write the function with the correct pattern and ensure no syntax errors.

Another possible issue: in the previous code, maybe the ketone pattern was incorrect. For example, using [CX3]=O without checking the neighboring atoms, which could include aldehydes or other groups.

But in the current approach, the pattern [CX3]([#6])([#6])=O ensures that the carbonyl carbon is connected to two other carbons, thus excluding aldehydes (which have only one carbon adjacent) and carboxylic acid derivatives.

Thus, the code should correctly identify diketones by counting two such groups.