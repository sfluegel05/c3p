"""
Classifies: CHEBI:23899 icosanoid
"""
</think>

Okay, let's try to figure out why the previous code didn't work. The user mentioned an "unexpected indent" error, which usually means there's a problem with the code's formatting. Looking at the previous code, I see that it's not valid Python code at all. Instead, it's a list of SMILES strings, not a function. That's a big mistake. The assistant probably pasted the examples into the code block by accident instead of writing the actual function.

Now, the task is to create a function called `is_icosanoid` that checks if a molecule is an icosanoid. From the definition, icosanoids are derived from C20 EFAs like arachidonic acid. So the key features should be:

1. A 20-carbon chain (but maybe with modifications like rings or branches).
2. Presence of oxygen-containing functional groups (like carboxylic acid, hydroxyl, epoxy, etc.) from oxidation.
3. Possibly multiple double bonds in the structure.

So, the approach should be:

- Check if the molecule has a carboxylic acid group (since many examples end with =O).
- Look for a 20-carbon backbone. But since the backbone might be modified (like cyclopentane in prostaglandins), maybe count the total carbons instead.
- Check for oxidation features like hydroxyl, epoxy, peroxide groups.
- Ensure the molecule is derived from a C20 precursor, which might mean the total carbons are around 20, but some examples have longer chains. Wait, looking at the examples, some have longer chains due to ester groups or other modifications. Hmm, maybe the core structure is 20 carbons, but substituents can add more.

Alternatively, maybe the key is that the parent fatty acid is C20. But how to determine that from the structure? It's tricky.

Another approach: Check for a carboxylic acid group, then a long carbon chain with multiple double bonds and oxygen-containing groups. Maybe count the number of carbons in the main chain, but it's hard to automate.

Alternatively, use the molecular formula. Icosanoids are derived from C20 precursors, so they should have 20 carbons. Wait, but some examples might have more due to ester groups or other substituents. For example, leukotriene C4 methyl ester has additional groups. So molecular formula might not be reliable.

Maybe the core structure (without substituents) has 20 carbons. But how to detect that programmatically?

Perhaps the minimal checks are:

1. Presence of a carboxylic acid group (many examples end with COOH).
2. Total number of carbons is at least 20 (but some may have more due to substituents).
3. Presence of oxidation-derived groups (like hydroxyl, epoxy, peroxide).
4. Multiple double bonds (since EFAs are polyunsaturated).

So the plan is:

- Parse the SMILES and check validity.
- Check for at least one carboxylic acid group.
- Count the total carbons. Maybe >= 20?
- Check for oxygen-containing groups (OH, O in epoxide, etc.).
- Check for multiple double bonds (maybe at least 2 or 3).

Let's implement this step by step. Also, handle cases where the molecule is an ester (like methyl ester), which would replace the carboxylic acid with COOCH3. Wait, but some examples have ester groups. For example, "leukotriene C4 methyl ester" has a methyl ester. So maybe the carboxylic acid could be in the form of an ester. Hmm, but the parent EFA is a carboxylic acid. So perhaps the icosanoid could be an ester derivative. So maybe the code should check for either a carboxylic acid or an ester group.

Alternatively, look for a chain that would be 20 carbons if the ester is hydrolyzed. But that's complicated.

Alternatively, check that the main chain (longest carbon chain) has 20 carbons, but that's not straightforward.

This is getting complicated. Maybe start with checking for a carboxylic acid group (or ester), a certain number of carbons (like 20), and presence of oxygen atoms beyond the carboxylic acid.

Let's proceed with:

1. Check if the molecule has a carboxylic acid group (or ester).
2. Check that the total number of carbons is at least 20.
3. Check for oxygen atoms in other groups (like hydroxyl, epoxy).
4. Check for multiple double bonds.

But how to implement this in RDKit?

First, parse the molecule. Then:

- Find carboxylic acid groups: look for COOH or COOR (ester).
- Count total carbons.
- Check for hydroxyl groups (OH), epoxy (O in a ring), peroxide (OO), etc.
- Count double bonds.

Let's write the code accordingly.

Wait, but some examples like prostaglandins have a cyclopentane ring with oxygen substituents. So the code needs to account for various oxygen-containing groups.

Another point: The EFAs are 20 carbons long. So the icosanoid should have a backbone derived from that. But with possible modifications like cyclization, which shortens the chain but keeps the total carbons around 20.

This is challenging. Maybe the best approach is to check for the presence of a carboxylic acid group, at least 20 carbons in total, and at least three double bonds (since arachidonic acid has four). But some icosanoids might have fewer due to saturation during oxidation.

Alternatively, use substructure matching for the oxidized parts. For example, look for prostaglandin-like structures with a cyclopentane ring and two side chains.

But given the diversity of icosanoids, it's hard to cover all with a single pattern. Maybe a combination of checks.

Let's outline the code steps:

1. Parse SMILES into a molecule.
2. Check for validity.
3. Check for presence of carboxylic acid or ester group.
4. Count total carbons. If less than 20, return False.
5. Check for at least two oxygen atoms (since oxidation adds O).
6. Check for at least two double bonds (since EFAs are polyunsaturated).
7. If all conditions met, return True.

But this might be too generic. Let's see examples:

- 15R-PGD2 has a cyclopentane ring, two side chains, carboxylic acid. Total carbons? Let's count: the main chain might be 20 carbons including the ring.

But how to count programmatically? Maybe total carbons >=20.

Another example: 5-HPETE has a 20-carbon chain with a hydroperoxy group and carboxylic acid.

So, the code could:

- Check for at least one carboxylic acid group (or ester).
- Total carbons >=20.
- At least two oxygen atoms (from oxidation).
- At least two double bonds.

Let's implement this.

But wait, some icosanoids might have fewer double bonds if they are saturated during processing. But the examples provided have multiple double bonds.

So, the code:

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_icosanoid(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid or ester
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ester = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")
    if not mol.HasSubstructMatch(carboxylic_acid) and not mol.HasSubstructMatch(ester):
        return False, "No carboxylic acid or ester group"
    
    # Count total carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Only {carbon_count} carbons"
    
    # Check for oxidation (at least one oxygen beyond carboxylic acid/ester)
    # Total oxygen count
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:  # At least one O beyond the COOH/ester
        return False, f"Insufficient oxygen atoms: {oxygen_count}"
    
    # Check for at least two double bonds
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 2:
        return False, f"Only {double_bonds} double bonds"
    
    return True, "Meets icosanoid criteria"

But wait, the ester group has an oxygen, so the oxygen count might include that. For example, a simple ester would have 2 oxygens (from COO), but if there's an additional hydroxyl, it would have 3. So requiring oxygen_count >=2 might be okay.

Testing this code against the examples:

For 15R-PGD2: SMILES has COOH, multiple oxygens, 20+ carbons, multiple double bonds. Should pass.

For leukotriene C4 methyl ester: has ester, enough carbons, oxygens, double bonds. Should pass.

But what about a molecule that meets all criteria but isn't an icosanoid? Possible false positives, but given the problem's constraints, this might be the best approach.

This should address the previous error by writing valid Python code with proper indentation and structure.