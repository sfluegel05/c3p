"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
</think>

Okay, I need to create a function called is_triterpenoid_saponin that checks if a given SMILES string represents a triterpenoid saponin. From the definition, a triterpenoid saponin is a terpene glycoside where the terpene part is a triterpenoid. 

First, I should break down the requirements. The molecule must have two main parts: a triterpenoid (which is a type of terpene with 30 carbons, derived from squalene) and a glycoside (sugar part attached via a glycosidic bond). 

So the steps I need to take are:

1. **Check for the presence of a triterpenoid skeleton.** Triterpenoids typically have a carbon skeleton with 30 carbons. They often have a complex structure with multiple rings, maybe like oleanane or ursane types. Detecting this via SMARTS might be tricky, but maybe I can look for a characteristic substructure or count the number of carbons.

2. **Identify glycosidic bonds.** Glycosides have sugar units attached via an oxygen atom (glycosidic bond). So I need to find at least one O-linked sugar group. The sugar part is usually a cyclic structure with multiple hydroxyl groups.

3. **Ensure the glycoside is attached to the triterpenoid.** The oxygen in the glycosidic bond should be connected to the triterpenoid part.

Possible approach using RDKit:

- Parse the SMILES into a molecule.
- Check if the molecule has a glycosidic group. A glycosidic bond is an oxygen atom connected to a sugar (which is a cyclic structure with multiple OH groups). Maybe use a SMARTS pattern for the glycosidic oxygen, like [OX2]-C1-C-O... but sugars can vary, so maybe look for a cyclic ether (like pyranose or furanose) with multiple hydroxyls.
- Check for the triterpenoid part. Since triterpenoids have 30 carbons, maybe check the molecular formula for C30HxOy... but some derivatives might have different H and O counts. Alternatively, look for common triterpenoid substructures, but that might be complex. Alternatively, check if the molecule has a large carbon skeleton (maybe by counting carbons and ensuring it's around 30, plus the sugar part).
- Ensure that the glycosidic oxygen is attached to the triterpenoid part.

Potential issues:

- Triterpenoid detection via SMARTS is difficult due to structural diversity. Maybe count the number of carbons in the aglycone (the part without the sugar). If the molecule has a sugar part, subtracting those carbons should leave around 30.
- The sugar part can vary in size (like glucose, xylose, etc.), so need to identify any sugar moiety.

Let me outline the steps in code:

1. Parse the SMILES into a molecule. Check validity.
2. Find all oxygen atoms connected to a sugar (maybe via a glycosidic bond). For this, perhaps look for oxygen atoms connected to a cyclic structure (like a hexose or pentose ring) with multiple hydroxyl groups.
3. Once a glycosidic oxygen is found, check if the other side of the bond is part of a triterpenoid structure.
4. For the triterpenoid part, check if the aglycone (molecule without the sugar parts) has around 30 carbons. Alternatively, check the total number of carbons in the molecule minus the sugar parts. But this might be complicated if there are multiple sugars.
5. Alternatively, check the molecular formula for a significant number of carbons (e.g., total carbons >=30 + sugar carbons). But this is not precise.

Alternative approach:

- Split the molecule into the aglycone and sugar parts. If the aglycone has 30 carbons (or close) and there's at least one sugar attached via O-glycosidic bond, then it's a triterpenoid saponin.

But splitting the molecule into aglycone and sugar parts programmatically is challenging. Maybe using the RDKit's Murcko decomposition or fragmenting around the glycosidic bond.

Another idea: Use the presence of a glycosidic bond (O connecting to a cyclic ether with multiple OHs) and check that the rest of the molecule (aglycone) has a triterpenoid structure. But how to define triterpenoid structure?

Alternative plan:

1. Check for the presence of at least one sugar group (cyclic structure with multiple OHs).
2. Check that the molecule has a large carbon skeleton (maybe total carbons >=30, but some saponins might have more if there are multiple sugars).
3. Verify that the sugar is attached via an O-glycosidic bond to the aglycone.

But the exact count might vary. For example, a triterpenoid (30 carbons) plus a hexose (6 carbons) would have 36 carbons. So total carbons >=30 might not be sufficient, but perhaps check that the majority of carbons are in the aglycone.

Alternatively, after removing the sugar parts, check if the remaining structure has around 30 carbons.

But how to remove the sugar parts? Maybe recursively remove any terminal sugar units.

This is getting complicated. Maybe start with:

- Check for the presence of at least one glycosidic bond (O linking to a sugar).
- Check that the aglycone (without sugars) has a molecular formula with 30 carbons.

But how to separate aglycone from sugars?

Alternatively, use a SMARTS pattern for the triterpenoid core. For example, oleanane type triterpenes have a specific fused ring system. But creating a SMARTS pattern for all possible triterpenoids is difficult.

Given the time constraints, perhaps a pragmatic approach:

- Check for the presence of at least one glycosidic group (sugar attached via O).
- Check that the molecule has a high number of carbons (e.g., >=30) which is typical for triterpenoids plus sugars.
- Additionally, check for common triterpenoid features, like multiple rings or a certain number of rings.

Alternatively, use the number of carbons in the aglycone. If the molecule has a glycosidic bond, subtract the sugar's carbons and see if the remaining is around 30.

But this requires identifying all sugar parts, which is non-trivial.

Perhaps the code can proceed as follows:

1. Check if the molecule has at least one glycosidic bond (O connected to a cyclic sugar).
2. Check that the molecular formula has at least 30 carbons (assuming that the aglycone is 30 and sugars add more).
3. Additionally, check for common triterpenoid features, like multiple rings or a certain substructure.

But this might not be accurate. For example, a diterpene (20 carbons) with a large sugar could have 30+ total carbons, leading to a false positive.

Alternative approach for triterpenoid detection: Triterpenoids are typically pentacyclic or tetracyclic. So check if the aglycone has at least four rings.

But again, separating the aglycone from the sugar is needed.

Another idea: Use the RDKit's Murcko decomposition to get the scaffold and check if it's complex enough (multiple rings) and has sufficient size.

But this is getting too vague. Given the examples provided, like Oleanolic acid derivatives, perhaps the presence of a characteristic triterpenoid skeleton can be detected via a SMARTS pattern.

For example, oleanolic acid has a specific pentacyclic structure. But creating a SMARTS pattern for that is challenging. Alternatively, look for a substructure with multiple fused rings.

Alternatively, check for a substructure that includes four or more fused rings, which is common in triterpenoids.

But without a specific SMARTS pattern, this is hard. Maybe count the number of rings in the molecule. Triterpenoids usually have 4 or 5 rings. But glycosides add the sugar rings, so the total might be higher. However, if the aglycone has 4-5 rings, that's a sign.

But again, separating aglycone rings from sugar rings is needed.

Given the time constraints, perhaps proceed with:

- Check for presence of glycosidic bonds (sugar moieties).
- Check that the molecule has a high number of carbons (>=30) and multiple rings.

But this is not foolproof.

Alternatively, use the following steps:

1. Check that the molecule has at least one sugar group (using a SMARTS for a cyclic ether with multiple OH groups).
2. Check that the molecule has a triterpenoid part, which can be approximated by a large carbon skeleton (e.g., >=30 carbons in the aglycone).

But how to compute the aglycone's carbons?

Maybe:

- Find all atoms in sugar groups (cyclic with multiple OHs).
- The aglycone is the rest of the molecule.
- Check if the aglycone has at least 30 carbons.

But implementing this requires identifying sugar groups, which is non-trivial.

Alternatively, use the following SMARTS for a glycosidic oxygen: [O;X2][C;!$(C=O)]@[C]. But this is not precise.

Alternatively, use a SMARTS pattern that matches a glycosidic bond, like oxygen connected to a carbon which is part of a ring (sugar).

Perhaps the code can proceed as:

- Use SMARTS to find oxygen atoms that are part of a glycosidic bond. For example, O connected to a carbon in a ring (sugar) and another carbon (aglycone).
- Once such oxygen is found, check the aglycone part (the part not in the sugar) for triterpenoid characteristics.

But how to split the molecule into aglycone and sugar parts.

Alternatively, for each glycosidic oxygen, split the molecule into two parts: the sugar and the rest. Then check if the rest has 30 carbons.

But splitting the molecule at the glycosidic bond would require breaking the bond and getting the fragments.

In RDKit, we can use the FragmentMol function to split into fragments. But this might split into multiple parts. For each glycosidic bond, split the molecule and check if one fragment is a sugar and the other is the aglycone.

But this is getting complicated. Let's try to outline code steps:

Function is_triterpenoid_saponin(smiles):

1. Parse the SMILES into a molecule. If invalid, return False.
2. Find all potential glycosidic bonds. These are oxygen atoms connected to a sugar (cyclic with multiple OHs) and another carbon (aglycone).
3. For each such bond, split the molecule into two parts: the sugar and the aglycone.
4. Check if the aglycone has around 30 carbons and a complex ring structure (e.g., at least 4 rings).
5. If any such split meets the criteria, return True.

But implementing this requires:

- Identifying glycosidic bonds.
- Splitting the molecule into aglycone and sugar.
- Checking the aglycone's properties.

Alternatively, use the following approach:

Check for the presence of at least one sugar group (using a SMARTS pattern for a cyclic sugar) and that the molecule has a triterpenoid part (large carbon skeleton).

But how to count carbons in the aglycone.

Alternatively, use the total number of carbons. Triterpenoid saponins usually have 30 (aglycone) + 1 or more sugars (each 5-6 carbons). So total carbons would be >=30 + 5 = 35. But this is a rough estimate.

But some examples may vary. For instance, the first example, dianversicoside E, has a SMILES that includes multiple sugar units. Let's count the carbons:

Looking at the SMILES, the triterpenoid part (aglycone) would have a large number of carbons, and the sugars add more. For example, if the aglycone is 30 carbons and there are two hexoses (12 carbons), total is 42. So checking for total carbons >=30 might not be sufficient, but maybe >=30 is a start.

But this approach could have false positives (e.g., a diterpene with multiple sugars).

Alternative idea: Use the RDKit's formula to get the number of carbons. If the molecule has >=30 C and has at least one glycosidic bond, return True.

But this is not accurate. For example, a molecule with 40 carbons but not a triterpenoid.

Given time constraints, perhaps proceed with:

- Check for presence of at least one glycosidic bond (sugar attached via O).
- Check that the molecule has a molecular weight typical of triterpenoid saponins (e.g., over 500 Da).
- Check that the number of carbons is at least 30.

But molecular weight can vary. Alternatively, use a combination of these factors.

Another approach: Look for the triterpenoid skeleton using a SMARTS pattern. For example, the oleanane skeleton can be represented with a specific pattern. But creating such a SMARTS is complex.

Alternatively, use a generic approach:

1. Check for glycosidic bonds (sugar moieties).
2. Check that the molecule has a large aglycone (>=30 carbons).

But again, how to split the aglycone.

Perhaps, in code:

- Find all the sugar rings (using SMARTS for cyclic ethers with multiple OHs).
- Remove those rings and count the carbons in the remaining structure.

But how to do that in RDKit.

Alternatively, for each atom in the molecule, determine if it's part of a sugar ring. If not, count those carbons as part of the aglycone.

But how to identify sugar rings. A possible SMARTS for a pyranose ring (6-membered with multiple OHs) might be:

[C]1[C][C][C][C][C]1 with multiple O attached.

But this is too vague. Alternatively, look for rings with multiple hydroxyl groups.

Alternatively, use the following SMARTS for a sugar-like ring:

[O;H1][C]1[C][C][C][C][C]1 (a 6-membered ring with an OH attached). But this is not comprehensive.

Given the time, perhaps proceed with the following code outline:

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_triterpenoid_saponin(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for at least one glycosidic bond (O connected to a sugar)
    # Sugar is a cyclic structure with multiple OHs
    # Use a SMARTS pattern for a glycosidic oxygen connected to a sugar ring
    glycosidic_o = Chem.MolFromSmarts('[O;X2][C;R]@[C]')
    if not mol.HasSubstructMatch(glycosidic_o):
        return False, "No glycosidic bond found"
    
    # Check molecular weight (triterpenoid saponins are usually >500)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} < 500)"
    
    # Check number of carbons in aglycone (approximate)
    # Assuming sugars are rings with 5-6 carbons, subtract from total
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # If total carbons >=30 (aglycone) + 6 (sugar) = 36, but this is a rough estimate
    if total_c < 30:
        return False, f"Total carbons ({total_c}) too low for triterpenoid"
    
    # Check for triterpenoid structure (complex ring system)
    # Count rings in the molecule, expecting multiple fused rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, f"Not enough rings ({n_rings}) for triterpenoid"
    
    return True, "Has glycosidic bond, sufficient size, and ring structure"

But this is a rough approach. The SMARTS pattern [O;X2][C;R]@[C] looks for an oxygen connected to a ring carbon, which might capture glycosidic bonds. The molecular weight check is arbitrary but might filter out small molecules. The total carbons check might be too lenient. The ring count includes both aglycone and sugar rings, which might not be accurate.

Testing with the example SMILES provided:

For example, Oleanolic acid 28-O-beta-D-glucopyranoside has a triterpenoid aglycone (30 carbons) plus a glucose (6 carbons), total 36. The code would check total_c >=30, which passes. The ring count would include the triterpenoid's rings (4-5) plus the sugar's ring, totaling >=5, so passes.

Another example, ginsenoside Rb2: the aglycone is a triterpenoid, plus multiple sugars. Total carbons would be higher than 30, passes. Ring count would be high.

Potential issues:

- Some molecules with multiple sugars but a non-triterpenoid structure might pass (false positives).
- The SMARTS pattern might miss some glycosidic bonds.
- The ring count might include sugar rings, leading to incorrect aglycone ring counts.

But given the time constraints and complexity, this approach might be a starting point. However, the previous code attempt had an indentation error. The user's last code attempt seems to have included a list of SMILES instead of code. So the actual error was a syntax error in the code structure.

In the current attempt, the code should avoid syntax errors. The code I proposed uses proper indentation and structure. Let's test it with the example SMILES.

For instance, the first example SMILES for dianversicoside E:

The code would parse the SMILES, find glycosidic bonds (since there are O-C links to sugar rings), check molecular weight (which is high), total carbons (>=30), and ring count (sufficient). So it would return True.

Thus, this code might work better than the previous attempt. However, there's a possibility of false negatives if the SMARTS pattern doesn't match the glycosidic bond correctly. Also, the ring count might not be sufficient. But given the problem's complexity, this is a reasonable approach.