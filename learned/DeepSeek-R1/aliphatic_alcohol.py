"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
</think>

Okay, let's see why the previous code had so many false positives and how to fix it. The main issue seems to be that the current approach only checks if the hydroxyl is attached to a non-aromatic carbon, but that's not enough. Aliphatic alcohols require the hydroxyl to be on a saturated (sp3) carbon in an aliphatic chain or ring. 

First, the SMARTS pattern [O;H1][C] matches any hydroxyl attached to any carbon, which includes things like carboxylic acids (where the OH is attached to a carbonyl carbon) or other functional groups. For example, in the false positive "rhamnetin 3-rhamninoside", there might be hydroxyls on carbons that are part of an aromatic system but the code mistakenly considers them aliphatic because the carbon isn't in a ring. Wait, no—the code checks if the carbon is aromatic. Hmm. Wait, maybe the problem is that the code only checks if the carbon is aromatic, but some hydroxyls are on non-aromatic carbons that are part of other functional groups, like esters or amides. For example, in "Arthpyrone D", there's a hydroxyl on a carbon adjacent to a carbonyl, making it part of an enol or similar structure, which isn't a typical aliphatic alcohol.

So the key is to ensure that the hydroxyl is attached to a carbon that's part of an aliphatic chain or ring and not adjacent to any electron-withdrawing groups (like carbonyls, nitros, etc.) that would make it a different functional group. Also, the carbon should be saturated (sp3 hybridized) to exclude enols or phenolic-like structures where the hydroxyl is conjugated.

To fix this, the SMARTS pattern should be adjusted to match hydroxyl groups attached to aliphatic carbons that are not part of any conjugated system. We can use a pattern like [C;!$(C=*)]-O, which matches a carbon not double-bonded to any atom. Also, check that the carbon has only single bonds and is not in an aromatic ring.

Wait, but how to ensure the carbon is aliphatic and not part of a conjugated system? Maybe using a SMARTS pattern that specifies a carbon with only single bonds and no adjacent electronegative atoms. Alternatively, check that the carbon with the hydroxyl is not adjacent to any double bonds or electronegative groups.

Another approach: The hydroxyl should be on a carbon that is not part of a carbonyl group (like carboxylic acid, ester, amide, etc.), not adjacent to a double bond (enol), and not in an aromatic ring. Also, the carbon should be sp3 hybridized.

So, modifying the SMARTS pattern to [C;X4;!$(C=O)]-[O;H1]. This would match a hydroxyl attached to a carbon with four bonds (sp3) and not part of a carbonyl group. However, RDKit's SMARTS might need to handle this carefully.

Alternatively, use two patterns: one to find hydroxyls, then check their environment. For each hydroxyl oxygen, check its neighboring carbon's hybridization and whether it's part of a carbonyl group.

Let me adjust the code:

1. Find all hydroxyl groups (O with one H and connected to a carbon).
2. For each such oxygen, get the attached carbon.
3. Check that the carbon is not aromatic (already in the code).
4. Check that the carbon is not part of a carbonyl group (C=O).
5. Check that the carbon's hybridization is sp3 (tetrahedral).

Implementing this would require checking each carbon for these conditions. Let's adjust the SMARTS pattern to exclude carbons adjacent to double bonds or electronegative atoms. Maybe:

Pattern for hydroxyl on aliphatic carbon: [C;!$(C=*);!$(C#*)]-[O;H1]

But this might not cover all cases. Alternatively, use a SMARTS that ensures the carbon is not part of a conjugated system.

Another approach: For each hydroxyl group, after finding the attached carbon, check that it has no double bonds and is not adjacent to any electronegative atoms (like O, N, S).

Alternatively, use the following steps:

- Check that the hydroxyl is attached to a carbon (C).
- The carbon must not be in an aromatic ring.
- The carbon must have at least one hydrogen (to exclude carbonyls where the C is double-bonded to O).
- The carbon must not be part of a carbonyl group (C=O) or other conjugated systems.

In RDKit, checking the hybridization might help. For example, if the carbon is sp3 hybridized, it's more likely to be an aliphatic alcohol.

So, modifying the code:

After finding the hydroxyl oxygen and the attached carbon, check:

1. Carbon is not aromatic (already in code).
2. Carbon's hybridization is SP3.
3. The carbon is not part of a carbonyl group (C=O).

Implementing this would reduce false positives where the hydroxyl is on a carbon adjacent to a carbonyl or in a conjugated system.

Let me adjust the code accordingly. Also, the previous code returns True on the first aliphatic hydroxyl found, which is correct, but the checks need to be stricter.

So, in the loop over matches:

for oxygen_idx, carbon_idx in matches:
    carbon = mol.GetAtomWithIdx(carbon_idx)
    if carbon.GetIsAromatic():
        continue
    # Check hybridization is SP3
    if carbon.GetHybridization() != Chem.HybridizationType.SP3:
        continue
    # Check that the carbon is not part of a carbonyl group
    # Look for adjacent C=O
    has_carbonyl = False
    for neighbor in carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 0:
            # Check if neighbor is double bonded to carbon
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                has_carbonyl = True
                break
    if has_carbonyl:
        continue
    # If all checks passed
    return True, "Hydroxyl group attached to aliphatic carbon"

This would exclude hydroxyls attached to carbons that are part of carbonyl groups (like carboxylic acids) and ensure the carbon is SP3.

Also, the initial SMARTS pattern might need to be adjusted to [O;H1][C;!$(C=*)] to avoid matching carbonyls.

Testing this should reduce false positives by ensuring the hydroxyl is on a saturated aliphatic carbon not adjacent to a carbonyl.

Another thing: The original code's SMARTS was [O;H1][C], which includes any carbon. By adding checks for SP3 and not being part of a carbonyl, we can filter out the false positives.

For example, in the false positive "Arthpyrone D", the hydroxyl might be attached to a carbon adjacent to a carbonyl, which this code would now exclude.

Similarly, for amino acids like "(2R,3R,4R,5S,6R)-2-Amino-4-hydroxy-3-methylpentanoic acid", the hydroxyl is on a carbon adjacent to a carbonyl (in the carboxylic acid group), so the code would skip it.

This should help reduce false positives while retaining true positives where the hydroxyl is on a simple aliphatic chain.