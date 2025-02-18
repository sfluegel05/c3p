"""
Classifies: CHEBI:26607 saturated fatty acid
"""
Now, write a new Python program.
</think>

Okay, I need to create a function called is_saturated_fatty_acid that checks if a given SMILES string represents a saturated fatty acid. Let's start by understanding what defines a saturated fatty acid. From the definition, it's a fatty acid with no carbon-carbon multiple bonds. Also, the examples provided all have a carboxylic acid group (O=O) attached to a carbon chain with only single bonds.

First, the function should check if the molecule has a carboxylic acid group. That's the -COOH group. So I'll look for the substructure where a carbon is double-bonded to an oxygen and single-bonded to another oxygen (like the SMARTS [CX3](=O)[OX2H1]).

Next, the molecule should have a carbon chain. Saturated means all single bonds between carbons. So I need to ensure there are no double or triple bonds in the carbon chain. But wait, the carboxylic acid group itself has a double bond (C=O), which is allowed. So I should exclude that from the check. So I need to find all the carbon-carbon bonds and make sure they are single.

But how do I separate the carboxylic acid part from the rest of the chain? Maybe after identifying the carboxylic acid group, check the rest of the molecule for any double or triple bonds. Alternatively, count all the bonds in the molecule except the ones in the carboxylic acid group and see if any are multiple.

Another approach: The entire carbon chain (excluding the carboxylic acid's C=O) should have only single bonds. So the molecule should have exactly one double bond (the one in the carboxylic acid). If there are more double bonds, it's unsaturated.

Wait, but some substituents might have other functional groups. For example, if there's a hydroxyl group (like in some examples), but the carbon chain itself is all single bonds. So the presence of other functional groups (like hydroxyl) doesn't affect the saturation of the fatty acid. The key is the carbon chain (the main backbone) being all single bonds.

But how to define the main chain? Maybe the longest carbon chain ending with the carboxylic acid. But that might complicate things. Alternatively, any carbon in the molecule (except the carboxylic acid's C=O) should not have double bonds.

So steps for the function:

1. Parse the SMILES into a molecule. If invalid, return False.
2. Check for the presence of a carboxylic acid group (-COOH). If not present, return False.
3. Check that there are no carbon-carbon double or triple bonds in the molecule, except for the one in the carboxylic acid group.
4. Ensure the carbon chain is sufficiently long. Fatty acids typically have at least 4 carbons (like butyric acid, CCCC(O)=O). But some examples have shorter chains, like propionic acid (CCC(O)=O) which has 3 carbons. Wait, looking at the examples: valproic acid is CCCC(CCC)C(O)=O, which might have branching. So the chain length might vary. However, the definition says "fatty acid", which generally refers to longer chains. But some of the examples include shorter ones like propionic (3 carbons), butyric (4), etc. So maybe the chain length is not a strict criterion here. However, perhaps the minimal is at least 2 carbons in the chain (since propionic has 3, but the carboxylic acid is two oxygens attached to a carbon). Wait, maybe the carboxylic acid is R-COOH, where R is an alkyl group. So the R group must be at least a methyl (but that would be acetic acid, which isn't in the examples). But the examples include propionic (3 carbons in total, R is ethyl). So maybe the function should not enforce a minimum chain length beyond having at least one carbon in the R group (so total of 2 carbons). However, some definitions of fatty acids require longer chains. But given the examples, perhaps we can't enforce a minimum length beyond having the carboxylic acid group. So maybe the chain length isn't part of the criteria here. The main check is the presence of the carboxylic acid and the absence of any other double bonds.

So, the key steps:

- Check for carboxylic acid group.
- Check that all other bonds (excluding the C=O in the acid) are single.

But how to implement that in RDKit?

First, find the carboxylic acid group. The SMARTS pattern for carboxylic acid could be [CX3](=O)[OX2H1]. But sometimes the hydrogen might be deprotonated (like in a salt), but the examples all have the protonated form. However, in SMILES, the acidic hydrogen might be written as part of the oxygen (like O), or as O-[Na+] for a sodium salt. But the problem states that the input is a SMILES string. So perhaps the function should check for the presence of the carboxylic acid group regardless of protonation state. Alternatively, look for the COOH group or its deprotonated form.

But maybe a better approach is to look for the presence of a carbon atom double-bonded to an oxygen (the carbonyl) and adjacent to an oxygen connected to a hydrogen (the hydroxyl). The SMARTS pattern for that could be [C;X3](=[O;X1])[O;X2][H]. But this might miss cases where the hydrogen is not explicitly present (like in a salt form). However, given the examples provided, perhaps assuming the protonated form is acceptable.

Alternatively, use a more general pattern: a carbon with a double bond to oxygen and another oxygen attached. For example, [CX3](=O)[OX2]. This would match both COOH and COO- (like in a salt). But then we need to ensure that the rest of the molecule doesn't have any other double bonds.

Wait, but the problem says "saturated fatty acid", which includes the COOH group. So the only double bond allowed is the one in the carboxylic acid. So the molecule should have exactly one double bond (the C=O in the acid), and no others.

So the steps would be:

1. Check that the molecule has exactly one double bond, which is part of the carboxylic acid group.
2. Check that there are no triple bonds.

But how to verify that the double bond is part of the carboxylic acid group?

Alternatively, check that all double bonds in the molecule are part of a carboxylic acid group.

So:

- Find all double bonds in the molecule.
- For each double bond, check if it's part of a carboxylic acid group.

If any double bond is not part of a carboxylic acid group, then it's not a saturated fatty acid.

Additionally, check that there are no triple bonds.

So code steps:

1. Parse the SMILES into a molecule. If invalid, return (False, "Invalid SMILES").
2. Check if the molecule has a carboxylic acid group. If not, return (False, "No carboxylic acid group").
3. Check all bonds in the molecule:
   a. Any triple bonds? If yes, return False.
   b. For each double bond, check if it's part of a carboxylic acid group. If any double bond is not part of such a group, return False.
4. If all checks pass, return True.

But how to check if a double bond is part of a carboxylic acid group?

The carboxylic acid group is -C(=O)OH. So the double bond is between C and O. So any double bond that is a C=O and adjacent to an O-H (or O connected to something else, like in a deprotonated form) would be part of the carboxylic acid.

Alternatively, the double bond in the carboxylic acid is between a carbon and an oxygen. So any double bond that is C=O is allowed (but only one such bond). Wait, but other functional groups like ketones (RCOR') would have C=O but are not carboxylic acids. So the presence of a carboxylic acid group is a separate check.

Wait, first, the molecule must have a carboxylic acid group. So step 2 ensures that. Then, any double bonds in the molecule must be part of that carboxylic acid group.

But how to find if there are any other double bonds.

So, after confirming that a carboxylic acid group exists, check that the total number of double bonds is exactly one (the one in the carboxylic acid group).

But wait, the carboxylic acid group has one double bond (C=O). So the molecule should have exactly one double bond. But if there are other double bonds (like in an unsaturated fatty acid), that's not allowed.

So code steps:

- After confirming the presence of a carboxylic acid group, check that the molecule has exactly one double bond (the C=O in the acid) and no triple bonds.

But wait, how to count the number of double bonds? Let's see:

mol = Chem.MolFromSmiles('CCCCCCCC(O)=O')  # octanoic acid
double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE]
len(double_bonds)  # should be 1 (the C=O)

Another example: a molecule with a double bond elsewhere, like C=CCCC(O)=O (unsaturated), would have two double bonds (the C=C and the C=O). So the check would be len(double_bonds) == 1.

But wait, in the carboxylic acid group, the C=O is one double bond. So if the molecule has any other double bonds (even C=O in another group), it's invalid.

So the code can check that the number of double bonds is exactly one, and that it's part of the carboxylic acid group.

But how to ensure that the double bond is part of the carboxylic acid group.

Alternatively, after finding the carboxylic acid group, check that all double bonds in the molecule are within that group.

So:

Find the atoms in the carboxylic acid group (C=O and adjacent O). Then, for each double bond in the molecule, check if it's between the C and O of the carboxylic acid.

But how to do that in RDKit.

Another approach: Use SMARTS to match the carboxylic acid group, and then check that any double bond in the molecule is part of that matched group.

For example, the SMARTS pattern for carboxylic acid could be [CX3](=O)[OX2H0,O-]. But perhaps a better approach is to match the entire carboxylic acid group and then check that all double bonds are part of that group.

But this could get complicated. Alternatively, after confirming that there's a carboxylic acid group, check that the molecule has exactly one double bond, which is part of that group.

But how to verify that the double bond is part of the carboxylic acid.

Alternatively, find all double bonds. For each, check if the atoms involved are C and O. If any double bond is between C and C, then it's invalid.

Wait, because in a saturated fatty acid, the only double bond is between C and O (the carboxylic acid). Any C=C or other double bonds (like C=N, etc.) would disqualify it.

So steps:

1. Check for presence of carboxylic acid group.
2. Check that all double bonds in the molecule are C=O.
3. Check that there are no triple bonds.

This would allow molecules with multiple C=O groups (like if there are two carboxylic acid groups), but according to the examples, the molecules have only one COOH group. But the definition of saturated fatty acid is about the absence of C-C multiple bonds, not the number of COOH groups.

Wait, the definition says "any fatty acid containing no carbon to carbon multiple bonds". So the key is that there are no C-C double or triple bonds. The presence of other double bonds (like C=O in other groups) is allowed, but according to the examples, all the examples have only the COOH group as the functional group. However, some examples have hydroxyl groups (like (20R)-20-hydroxyhenicosanoic acid), which have an -OH group but no double bonds in the carbon chain.

So the main check is that there are no C-C double or triple bonds. The presence of other functional groups (like hydroxyl, ester, etc.) is allowed as long as they don't introduce C-C multiple bonds.

Wait, but the definition says "no carbon to carbon multiple bonds". So any double or triple bond between two carbon atoms is forbidden. The C=O in the carboxylic acid is allowed.

So the correct approach is:

1. Check that the molecule has a carboxylic acid group.
2. Check that there are no carbon-carbon double or triple bonds in the molecule.

So, code steps:

- Parse SMILES.
- Check for carboxylic acid group (using SMARTS pattern).
- Iterate through all bonds in the molecule. For each bond:
   - If it's a double or triple bond between two carbons, return False.

But how to implement that.

So, in RDKit, for each bond in mol.GetBonds():

bond_type = bond.GetBondType()
if bond_type in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
   begin_atom = bond.GetBeginAtom()
   end_atom = bond.GetEndAtom()
   if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
       return False, "Contains carbon-carbon multiple bond"

But wait, this would check for any C-C double or triple bonds, which is exactly what we need.

But also, the carboxylic acid group has a C=O double bond, which is allowed. So the code should not flag that.

So, the code would:

- Check for presence of carboxylic acid group.
- Check that there are no C-C double or triple bonds.

So the code would be:

def is_saturated_fatty_acid(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (-COOH or COO-)
    # SMARTS pattern for carboxylic acid: [CX3](=O)[OX2H1]
    # But RDKit might not handle H properly, so use [CX3](=O)[OX2H0,O-]?
    # Alternative: check for a carbon with a double bond to O and adjacent to an O
    carb_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carb_acid_pattern):
        # Check for deprotonated form (COO-)
        carb_acid_pattern_deprotonated = Chem.MolFromSmarts("[CX3](=O)[O-]")
        if not mol.HasSubstructMatch(carb_acid_pattern_deprotonated):
            return False, "No carboxylic acid group"
    
    # Check for any carbon-carbon double or triple bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
            begin = bond.GetBeginAtom()
            end = bond.GetEndAtom()
            if begin.GetAtomicNum() == 6 and end.GetAtomicNum() == 6:
                return False, "Contains carbon-carbon multiple bond"
    
    return True, "Saturated fatty acid with no carbon-carbon multiple bonds"

Wait, but in the carboxylic acid group, the C=O is a double bond between C and O. So that's allowed. The code checks for bonds between two carbons. So that's correct.

But what about other double bonds, like between C and N or O? For example, if there's a ketone group (C=O) elsewhere, but that's not part of the carboxylic acid. However, according to the definition, as long as there are no C-C multiple bonds, it's allowed. But the presence of other functional groups (like ketones) would make it not a fatty acid. Wait, but the definition says "any fatty acid". So maybe the code should also check that the only functional group is the carboxylic acid. But the examples include molecules with hydroxyl groups (like 20-hydroxyhenicosanoic acid). So other functional groups are allowed as long as there are no C-C multiple bonds.

Thus, the code should only check for the presence of a carboxylic acid group and the absence of C-C multiple bonds.

But wait, the code as written would return True for a molecule with a carboxylic acid group and another C=O group (like a ketone), as long as there are no C-C double bonds. For example, if the SMILES is CC(=O)CCCC(O)=O (a molecule with a ketone and a carboxylic acid group), the code would allow it. But is such a molecule considered a fatty acid? Probably not, as fatty acids are typically straight-chain carboxylic acids. However, according to the problem's definition, the key is no C-C multiple bonds. So the code would incorrectly classify such a molecule as a saturated fatty acid. But the examples provided don't include such cases, so maybe the code is sufficient for the given examples.

But perhaps the definition of a fatty acid requires that the molecule is a carboxylic acid with a carbon chain (aliphatic), possibly with substituents but no other functional groups. However, the examples include substituents like hydroxyl and methyl groups. So the code should not check for the absence of other functional groups, only the presence of the carboxylic acid and the absence of C-C multiple bonds.

So the code as written would correctly identify the examples given.

Testing with octanoic acid (CCCCCCCC(O)=O):

- The SMILES has a carboxylic acid group (matches the pattern).
- All bonds are single except the C=O in the acid. So the code returns True.

Testing a molecule with a C=C bond (like oleic acid, which is unsaturated), the code would detect the C=C and return False.

Another test case: 4-oxopentanoic acid (SMILES: CC(=O)CCC(O)=O). This has two C=O groups (one ketone, one carboxylic acid). The code would check for C-C bonds. The ketone's C=O is between C and O, so no C-C bonds. So the code would return True. But according to the examples, 4-oxopentanoic acid is included as a valid example. Wait, looking back at the examples, yes, 4-oxopentanoic acid is listed as an example. So this molecule has a ketone group (C=O) in addition to the carboxylic acid. But according to the code, since there are no C-C double bonds, it's considered a saturated fatty acid. Which aligns with the problem's examples. So the code would correctly classify it.

So the code seems to handle that case.

Now, the code's first part checks for the presence of a carboxylic acid group. The SMARTS pattern [CX3](=O)[OX2H1] matches the protonated form. The deprotonated form (COO-) would be matched by [CX3](=O)[O-]. However, in the examples, the SMILES strings are written with the protonated form (O)=O. But if the input SMILES is deprotonated (like CCCC(=O)[O-]), the code would check for that.

But perhaps the code should use a more general SMARTS pattern that matches both forms. Let's adjust the code to check for either the protonated or deprotonated form.

Alternatively, use a SMARTS pattern that matches the carboxylic acid group regardless of protonation. The SMARTS pattern [CX3](=O)[OX2H1,O-] would match both COOH and COO-.

But in RDKit, the H in the SMARTS might not be handled correctly if the molecule is not properly protonated. So perhaps using [CX3](=O)[OX2H0,O-] would match both.

Alternatively, use two separate checks: first for the protonated form, then for the deprotonated.

In the code above, the first check is for [CX3](=O)[OX2H1], and if that's not found, check for [CX3](=O)[O-].

But maybe there's a better way. Let's use a single SMARTS pattern that matches either.

The SMARTS pattern for carboxylic acid (protonated or deprotonated) can be written as [CX3](=O)[OX2H1,O-].

But in RDKit, the H1 might not be present if the molecule is not protonated. For example, in the SMILES "CC(=O)O", the oxygen has two single bonds (to C and H). But if the molecule is written as "CC(=O)[O-]", then the oxygen has a negative charge.

So the code can check for either case.

But how to write that in SMARTS.

In RDKit's SMARTS, [OX2H1] matches an oxygen with two bonds and one hydrogen. [O-] matches an oxygen with a negative charge. So the SMARTS pattern can be "[CX3](=O)[OX2H1,O-]".

But how to combine those in a single pattern.

Yes, the SMARTS would be "[CX3](=O)[OX2H1,O-]".

So in the code:

carb_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
if not mol.HasSubstructMatch(carb_acid_pattern):
    return False, "No carboxylic acid group"

This would match both protonated and deprotonated forms.

But perhaps the code should also check that the oxygen is connected to the carbonyl carbon. Wait, the SMARTS [CX3](=O)[OX2H1,O-] would match a carbonyl carbon (C=O) with an adjacent oxygen (either O-H or O-).

Yes, that's correct.

So, the code should first check for the presence of that pattern.

Then, check for any C-C double or triple bonds.

Thus, the code would be:

from rdkit import Chem
from rdkit.Chem import AllChem

def is_saturated_fatty_acid(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (protonated or deprotonated)
    carb_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    if not mol.HasSubstructMatch(carb_acid_pattern):
        return False, "No carboxylic acid group"
    
    # Check for any carbon-carbon double or triple bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                return False, "Contains carbon-carbon multiple bond"
    
    return True, "Saturated fatty acid with no carbon-carbon multiple bonds"

Testing this code with the example SMILES:

For octanoic acid (CCCCCCCC(O)=O):

- The SMILES is parsed into a molecule.
- The carboxylic acid pattern matches (C=O and adjacent O with H).
- All bonds are single except the C=O. So no C-C double bonds. Returns True.

For 4-oxopentanoic acid (CC(=O)CCC(O)=O):

- The carboxylic acid group is present (O=O).
- The molecule has two double bonds (C=O in the ketone and the carboxylic acid). But the code checks for C-C double bonds. The ketone's C=O is C-O, so no C-C bonds. So returns True. Which is correct according to the example.

For a molecule like C=CCCC(O)=O (unsaturated), the code would find the C=C bond (C-C double bond), returns False.

Another example: (20R)-20-hydroxyhenicosanoic acid (C[C@@H](O)CCCCCCCCCCCCCCCCCCC(O)=O). The code checks for carboxylic acid and no C-C multiple bonds. The hydroxyl group is allowed. Returns True.

Thus, this code should correctly classify the examples provided.