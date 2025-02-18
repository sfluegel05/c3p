"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
</think>

Okay, let's try to figure out why the previous attempt failed. The error message says "unexpected indent" in line 1. Looking at the previous code, it seems like the user pasted a SMILES string that was part of the examples as the previous code, which is invalid Python. That's a formatting mistake. The actual code probably wasn't properly enclosed in triple backticks, leading to the SMILES examples being interpreted as code.

Now, focusing on the task: classifying glycosaminoglycans. The definition is "Any polysaccharide containing a substantial proportion of aminomonosaccharide residues." So, the key features are:

1. It's a polysaccharide (long chain of monosaccharide units).
2. Contains a significant number of amino sugars (like glucosamine or galactosamine).

Steps to approach:
- Check if the molecule is a polysaccharide. This can be tricky, but perhaps using the presence of multiple saccharide units connected by glycosidic bonds.
- Identify aminomonosaccharide residues. These are sugars where a hydroxyl group is replaced by an amino group. The SMARTS pattern for an amino sugar could be something like [C;R][NH2] in a ring structure.
- Ensure that the proportion of amino sugars is "substantial," which is vague. Maybe check that at least 20% of the saccharide units have amino groups.

Possible RDKit methods:
- Use SMARTS to find glycosidic bonds (O linking two carbohydrate rings).
- Use SMARTS to find amino groups attached to the saccharide rings.
- Count the ratio of amino-containing saccharides to total saccharides.

But wait, RDKit might not directly identify all saccharide units. Alternative approach: find all cyclic ethers (like pyranose rings) and check for adjacent amino groups.

Another challenge: SMILES strings for polysaccharides can be complex with multiple repeating units. The examples provided seem to have complex structures with multiple rings and substituents, possibly including amino groups.

Potential code steps:
1. Parse the SMILES into a molecule.
2. Find all rings (assuming saccharide units are cyclic).
3. For each ring, check if there's an amino group attached (NH2 or NHR where R could be acetyl, etc.).
4. Count the number of amino-containing rings vs. total rings.
5. If the ratio is above a threshold (e.g., 20%) and there are multiple rings (polysaccharide), return True.

But this might not be accurate. Let's try writing the code with these steps.

Wait, the examples given don't look like typical glycosaminoglycans (e.g., hyaluronic acid, heparin). The provided SMILES have complex structures with many non-carbohydrate parts. Maybe the examples are incorrect, or perhaps the approach needs to be adjusted.

Alternative idea: Look for the presence of both carbohydrate features (multiple oxygen atoms in rings, glycosidic bonds) and amino groups attached to sugar rings.

So, the code could:
- Check for multiple glycosidic bonds (O connecting two carbons in rings).
- Check for amino groups (NH2, NHC(=O)R, etc.) attached to the carbohydrate rings.
- Ensure that a significant number of the rings have such amino groups.

Implementing this with RDKit's SMARTS patterns:

Glycosidic bond pattern: something like C-O-C in a chain connecting two rings.

Amino sugar pattern: a ring (like six-membered) with an NH2 or N-acetyl group attached.

But creating accurate SMARTS for this is challenging. Let's draft possible patterns.

First, find all rings in the molecule. Then, for each ring, check if any atom in the ring is a carbon connected to an amino group.

Alternatively, use a SMARTS pattern like:

[A]a1aaaaa1 where [A] is a carbon, and one of the substituents is an amino group.

But maybe better to look for a six-membered ring (like glucosamine) with an NH2 or N-acetyl group.

For example, the SMARTS for glucosamine could be:

[C]1[C][C][C][C]([NH2])[C]O1

But this might be too specific. Alternatively, look for any ring (size 5 or 6) with an amino group attached.

So, the code steps would be:

1. Parse the molecule.
2. Find all rings.
3. For each ring, check if any atom in the ring is a carbon connected to an amino group (NH2, NHCOCH3, etc.).
4. Count the number of such rings (amino_sugar_rings).
5. Total rings (total_sugar_rings) – but need to ensure they are sugar rings (maybe by checking for multiple oxygen atoms in the ring? Or adjacent hydroxyl groups?).
6. If amino_sugar_rings / total_sugar_rings >= 0.2 and total_sugar_rings >= 2 (indicating a polysaccharide), return True.

This is a simplified approach but might work for the given examples.

Potential code:

Use RDKit's ring detection, then for each ring, check for adjacent amino groups. Also, check that there are multiple rings (polysaccharide).

But how to distinguish sugar rings from other rings? Maybe by the presence of multiple oxygen atoms in the ring (like in pyranose).

Alternatively, use a SMARTS pattern for a sugar ring with an amino group.

But this is getting complicated. Let's proceed step by step.

First, the code structure:

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Find all rings
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) < 2:
        return False, "Not enough rings for a polysaccharide"
    
    amino_sugar_count = 0
    # Pattern for amino group attached to a ring carbon
    amino_pattern = Chem.MolFromSmarts('[NH2,nH1;!$(NC=O)]')  # Adjust as needed
    for ring in rings:
        # Check if the ring is a sugar (maybe has oxygen atoms in the ring)
        oxygen_in_ring = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
        if not oxygen_in_ring:
            continue  # Skip non-sugar rings
        # Check if any atom in the ring has an amino group
        for a in ring:
            atom = mol.GetAtomWithIdx(a)
            if atom.GetAtomicNum() == 6:  # Carbon
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 7:  # Nitrogen
                        # Check if it's an amino group (primary or secondary)
                        if any(bond.GetBondType() == Chem.BondType.SINGLE for bond in neighbor.GetBonds()):
                            amino_sugar_count +=1
                            break  # One per ring
                if amino_sugar_count >0:
                    break  # Move to next ring
    
    if amino_sugar_count ==0:
        return False, "No amino sugars detected"
    
    # Check proportion
    total_sugar_rings = len([ring for ring in rings if any(mol.GetAtomWithIdx(a).GetAtomicNum() ==8 for a in ring)])
    if total_sugar_rings ==0:
        return False, "No sugar rings found"
    
    ratio = amino_sugar_count / total_sugar_rings
    if ratio >= 0.2:
        return True, f"Amino sugar ratio {ratio:.2f} meets threshold"
    else:
        return False, f"Amino sugar ratio {ratio:.2f} too low"

But this code has several assumptions. For example, it counts any ring with an oxygen as a sugar ring, which might not be accurate. Also, the amino pattern detection is simplistic.

Testing with the example SMILES provided (e.g., Binchamide A), which may not be a glycosaminoglycan. Wait, the examples given might not be correct. For instance, Binchamide A's SMILES doesn't look like a polysaccharide. So perhaps the examples are misleading, making the problem harder.

Given the complexity, the code might return (None, None) if the task is too hard, but the user wants an attempt.

Another approach: Use the presence of multiple glycosidic linkages (O connecting two carbons in different rings) and amino groups in some of those rings.

But detecting glycosidic bonds accurately requires specific patterns.

Alternatively, use the number of carbohydrate monomers. RDKit has a module for carbohydrate handling, but it's not straightforward.

Given time constraints, perhaps proceed with the code above, adjusting the patterns and thresholds.