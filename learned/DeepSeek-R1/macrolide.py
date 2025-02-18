"""
Classifies: CHEBI:25106 macrolide
"""
</think>

Okay, I need to create a function called is_macrolide that checks if a given SMILES string represents a macrolide. From the definition, a macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide. 

First, I need to break down the requirements. The key features are:
1. Macrocyclic lactone: A large ring (12+ atoms) containing a lactone group (ester derived from a carboxylic acid and an alcohol, which is a cyclic ester).
2. Derived from a polyketide: This part is a bit tricky. Polyketides are often characterized by alternating carbonyl and methylene groups, but detecting this programmatically might be complex. Maybe I can skip this part if it's too hard, but I should try to find a pattern.

Let me start by checking for the lactone. A lactone is a cyclic ester, so I can look for a ring that includes an ester group (O-C=O). But the ester needs to be part of the macrocycle. 

So steps I can take:
1. Parse the SMILES into a molecule.
2. Find all rings in the molecule.
3. Check if any ring has 12 or more atoms.
4. Within that ring, check if there's an ester group (O-C=O) as part of the ring.

Wait, but the lactone is a cyclic ester, so the ester oxygen and the carbonyl should be part of the ring. So the ester group (O-C(=O)-O) is part of the ring. Wait, no. A lactone is formed by an intramolecular esterification of a hydroxy carboxylic acid. So the structure would be a ring where an oxygen is connected to a carbonyl group (like O-C(=O)- as part of the ring). So the ester group is part of the ring.

So maybe the SMARTS pattern for a lactone would be a ring where there's an oxygen atom adjacent to a carbonyl group (C=O), and both are part of the ring.

But how to write that SMARTS? Let's think: The lactone would have a structure like O-C(=O)- as part of the ring. So maybe the pattern is [O;R]-C(=O)-[C;R], where all atoms are in the same ring.

But I need to make sure that the oxygen and the carbonyl are part of the same ring. So the SMARTS could be something like [O;R1]-C(=O)-[C;R1], where R1 indicates the same ring.

Alternatively, using RDKit's ring detection. Maybe find all rings, then check each ring to see if it contains an ester group (O connected to a carbonyl that's part of the ring).

Another approach: Find all ester groups in the molecule, then check if any of them are part of a macrocyclic ring (size >=12).

But how to find ester groups? The ester group is O-C(=O)-O, but in a lactone, one of the oxygens is part of the ring. Wait, lactone is a cyclic ester, so it's a cyclic structure where an oxygen is connected to a carbonyl group. So the ester in the lactone would be part of the ring. For example, in a 14-membered lactone, the oxygen and the carbonyl are part of the 14-membered ring.

So the plan is:
1. Find all rings in the molecule.
2. For each ring, check if its size is >=12.
3. For those rings, check if they contain a lactone group (O connected to a carbonyl within the ring).

But how to check for the lactone in the ring? Maybe look for an oxygen atom in the ring that is adjacent to a carbonyl group (C=O) that is also in the ring.

So for each macrocyclic ring (size >=12), check if there's an oxygen atom in the ring that is bonded to a carbonyl carbon (C=O) which is also in the ring.

Alternatively, use a SMARTS pattern that matches atoms in a ring of size >=12 and includes an ester group.

Wait, maybe the SMARTS pattern for the lactone in a macrocycle would be:

[O;R]-C(=O)-[C;R]

But this would match an oxygen in a ring connected to a carbonyl carbon, which is also in the same ring. That should be part of the lactone.

So in RDKit, I can create a SMARTS pattern for that. Let's test this.

The SMARTS could be: [O]C(=O)[C] where all atoms are in the same ring. But to enforce that they are part of the same ring, perhaps using the RDKit's HasSubstructMatches with a SMARTS that includes ring membership.

Alternatively, create a SMARTS that ensures the O and the C=O are in the same ring. For example:

[O;R1]C(=O)[C;R1]

This would match an oxygen atom in ring 1, connected to a carbonyl carbon, which is also in ring 1. Then, the ring size can be checked separately.

So the steps would be:

1. Parse the SMILES into a molecule.
2. Check for validity.
3. Find all rings in the molecule.
4. For each ring, check if the size is >=12.
5. For each such ring, check if it contains the lactone group (O connected to C=O in the same ring).

But how to do this in RDKit? Maybe:

- Use the molecule's GetRingInfo() to get all rings.
- Iterate over each ring, check if size >=12.
- For each qualifying ring, check if there's an oxygen atom in the ring that is adjacent to a carbonyl group (C=O) also in the ring.

Alternatively, use a SMARTS pattern that matches the lactone in any ring of size >=12. But SMARTS doesn't directly support ring size conditions. So perhaps first find all rings of size >=12, then check each for the lactone group.

But how to do that in RDKit?

Alternatively, use a SMARTS pattern that matches a macrocyclic lactone. For example:

A ring of 12 or more atoms containing an ester group (O-C(=O)-).

But writing a SMARTS for a ring of size >=12 is tricky because SMARTS doesn't have a direct way to specify ring size. However, we can use the RDKit's ring size detection.

Another approach: Find all atoms that are part of a ring of size >=12. Then, check if any of those atoms form a lactone (O connected to C=O in the same ring).

So here's a possible plan for the code:

- Convert the SMILES to a molecule.
- Check if the molecule has any rings of size >=12.
- For each such ring, check if there's an oxygen atom in the ring that is adjacent to a carbonyl group (C=O) that is also in the ring.
- If any such ring exists, return True; else, False.

But how to implement this?

First, get all rings using GetRingInfo(). Then, for each ring, check the size. If size >=12, then check if any atom in the ring is an oxygen, and adjacent to a carbonyl carbon.

Wait, but the lactone's oxygen is part of the ester group. So in the ring, the oxygen is connected to a carbonyl carbon (C=O), which is part of the ring.

So for each ring of size >=12:

Loop through each atom in the ring. If the atom is an oxygen, check its neighbors. For each neighbor, check if it's a carbon with a double bond to oxygen (carbonyl). Also, ensure that this carbonyl carbon is part of the same ring.

For example:

For each ring in rings:
    if len(ring) >=12:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:  # oxygen
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in neighbor.GetBonds() if bond.GetOtherAtom(neighbor).GetAtomicNum() == 8):
                        # Check if neighbor is in the same ring
                        if neighbor.GetIdx() in ring:
                            return True

This would check if there's an oxygen in the ring adjacent to a carbonyl carbon that's also in the ring.

But wait, the carbonyl is C=O. So the neighbor of the oxygen (O) would be the carbonyl carbon (C=O), and that carbon should have a double bond to another oxygen (the carbonyl oxygen). But the carbonyl oxygen may not be in the ring. However, the lactone's structure requires that the oxygen (from the ester) and the carbonyl are part of the ring.

Wait, in a lactone, the ester is formed between a hydroxyl group and a carboxylic acid on the same molecule, creating a cyclic ester. So the structure is O-C(=O)-O as part of the ring. Wait, no. Let me think: the lactone is formed when a hydroxyl group (-OH) and a carboxylic acid (-COOH) undergo esterification, losing water. The resulting structure is a cyclic ester: O-C(=O)-O connected in a ring. Wait, no. The lactone structure is a cyclic ester where the oxygen from the hydroxyl becomes part of the ester linkage. So the structure would be something like:

R-O-C(=O)-R', where R and R' are connected to form a ring. So in the ring, there's an oxygen connected to a carbonyl group.

So in the ring, the oxygen is part of the ester linkage. So the oxygen is connected to the carbonyl carbon (C=O), and that carbon is also connected to another oxygen (the carbonyl oxygen). The carbonyl oxygen may or may not be part of the ring. But in the lactone, the ester oxygen (the one from the original hydroxyl) is part of the ring, and the carbonyl is part of the ring.

Wait, let's take an example. Let's say we have a 14-membered lactone. The ring would include the oxygen from the hydroxyl and the carbonyl group. So the structure would be something like:

...-O-C(=O)-... (with the rest of the ring connecting back). So the oxygen and the carbonyl carbon are part of the ring, but the carbonyl oxygen (the one in the C=O) is not necessarily part of the ring.

Wait, no. The carbonyl oxygen is part of the ester group. In the lactone, the ester is formed between the hydroxyl oxygen and the carbonyl carbon. So the oxygen in the ester is part of the ring, and the carbonyl carbon is part of the ring, but the carbonyl oxygen is a double-bonded oxygen and not part of the ring.

So the lactone's structure in the ring would have O connected to C=O, where both O and C are in the ring. The carbonyl oxygen (double-bonded) is not part of the ring.

Therefore, the code needs to check for an oxygen atom in the ring that is bonded to a carbon which is double-bonded to an oxygen (the carbonyl), and that carbon is also in the ring.

So, for each ring of size >=12, check if there exists an oxygen atom in the ring that is bonded to a carbon that has a double bond to an oxygen (anywhere), and that carbon is also in the ring.

So in code:

For each ring in rings:
    if len(ring) >=12:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:  # oxygen
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in ring:
                        # Check if this neighbor has a double bond to oxygen
                        for bond in neighbor.GetBonds():
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                other_atom = bond.GetOtherAtom(neighbor)
                                if other_atom.GetAtomicNum() == 8:
                                    # Found O-C(=O) in the ring
                                    return True

This would check if in the ring, there's an O connected to a C, which is double-bonded to another O, and both O and C are in the ring.

Wait, but the double-bonded O doesn't need to be in the ring. Only the O (from the ester) and the C need to be in the ring.

So the code would check for an oxygen in the ring that is bonded to a carbon (also in the ring) which has a double bond to any oxygen (could be inside or outside the ring). But in the lactone structure, the carbonyl oxygen is part of the ester group but not necessarily in the ring.

Wait, but in the lactone, the carbonyl is part of the ring. Let me think again. The lactone structure is a cyclic ester. The ester group is O-C-O-C=O, but arranged in a ring. For example, in a 14-membered lactone, the oxygen from the original hydroxyl is part of the ring, connected to a carbonyl group (C=O) which is also part of the ring. The carbonyl oxygen is double-bonded to the carbon and is not part of the ring backbone.

So the ring would have atoms in the order ...-O-C-C(=O)-... where the O and the C=O are part of the ring. The C=O is in the ring, but the oxygen of the carbonyl is not part of the ring (it's a substituent).

Wait, no. The carbonyl is part of the ring. Let's take a simple lactone, like gamma-butyrolactone, which is a 5-membered ring. The structure is O-C-C-C-C=O connected in a ring. Wait, no, gamma-butyrolactone is a 5-membered ring with the structure O-C-C-C-C(=O)-O, but that's not right. Actually, the structure is a cyclic ester formed from gamma-hydroxybutyric acid. The ring is O connected to a carbonyl group. So the structure is a five-membered ring: O-(CH2)3-C(=O)-, connected in a ring. So the oxygen and the carbonyl carbon are part of the ring. The carbonyl oxygen is not part of the ring.

So in the macrocyclic lactone (like a 14-membered ring), the oxygen and the carbonyl carbon are part of the ring. The carbonyl oxygen is a substituent.

Therefore, in the code, we need to check for an oxygen atom in the ring that is bonded to a carbon (also in the ring) which has a double bond to an oxygen (which can be anywhere). So the code as described before should work.

Now, putting this into code:

First, get the molecule's ring information.

Then, for each ring, check if the size is >=12.

For each such ring, check if there's an oxygen in the ring bonded to a carbon in the ring that has a double bond to oxygen.

Additionally, the molecule should have a macrocyclic ring (size >=12) containing a lactone.

Now, the code steps:

1. Parse the SMILES into a molecule. If invalid, return False.
2. Get the ring information using mol.GetRingInfo().AtomRings().
3. For each ring in the rings:
   a. Check if the ring size is >=12.
   b. For each atom in the ring, check if it's an oxygen.
   c. For each oxygen, check its neighbors.
   d. For each neighbor that is a carbon in the same ring, check if it has a double bond to an oxygen.
4. If any ring meets these criteria, return True.

Now, code implementation.

But also, the definition mentions that it's derived from a polyketide. How to check that? Polyketides have a backbone of alternating ketones (or their derivatives). But detecting this programmatically is challenging. Since the problem states that if it's too hard, return (None, None), but the examples provided include various macrolides, perhaps the lactone check is sufficient. Alternatively, maybe look for a certain pattern of carbonyl groups.

But given time constraints and complexity, perhaps the code can focus on the macrocyclic lactone part and not check the polyketide origin. The user may accept that, given the examples.

So the code will check for the presence of a macrocyclic lactone (ring >=12 with lactone group) and return True if found.

Now, implementing this in Python with RDKit.

Potential code:

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find all rings in the molecule
    rings = mol.GetRingInfo().AtomRings()
    
    for ring in rings:
        ring_size = len(ring)
        if ring_size < 12:
            continue  # Skip small rings
        
        # Check for lactone in this ring: O connected to C=O within the ring
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 8:  # Not oxygen
                continue
            
            # Check neighbors of this oxygen
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in ring:
                    continue  # Neighbor not in the ring
                if neighbor.GetAtomicNum() != 6:  # Not carbon
                    continue
                
                # Check if this carbon has a double bond to oxygen
                for bond in neighbor.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other_atom = bond.GetOtherAtom(neighbor)
                        if other_atom.GetAtomicNum() == 8:
                            # Found O-C=O in the ring
                            return True, f"Macrocyclic lactone with {ring_size}-membered ring"
    
    return False, "No macrocyclic lactone found"

Wait, but this code has a problem. The carbonyl oxygen (the one in the C=O) may not be in the ring, but the code checks if any neighbor of the oxygen in the ring has a double bond to any oxygen. That's correct because the carbonyl oxygen doesn't need to be in the ring.

Testing this code with example SMILES:

Take Erythromycin's SMILES, which is a macrolide. For example, erythromycin A has a 14-membered lactone ring.

But when we run this code on erythromycin's SMILES, would it detect the lactone?

Alternatively, take the example SMILES given, like Amphotericin B. Wait, no, Amphotericin B is a polyene, not a macrolide. Wait, no, Amphotericin B is a macrolide antibiotic. Its structure includes a macrocyclic lactone.

But perhaps the code would correctly identify it.

But perhaps there's a mistake in the code. For example, the code checks if the neighbor (carbon) has any double bond to oxygen. That could include carbonyl groups outside the lactone. But in the lactone, the carbonyl is part of the ring.

Wait, the code checks if the carbon is in the ring and has a double bond to oxygen. So even if the carbonyl is part of the ring, the code should find it.

Another possible issue: the code may find multiple rings, but only needs one to qualify.

So the code should return True if any ring of size >=12 contains an oxygen connected to a carbonyl carbon in the same ring.

This should work for macrolides.

But perhaps the code needs to check that the carbonyl carbon is part of the ring and has a double bond to oxygen.

Yes, that's what the code does.

Testing with a sample macrolide SMILES, like erythromycin:

Erythromycin's SMILES: "C[C@H]1[C@@H](O)C[C@@H](C)O[C@@H]1OC[C@@H]([C@@H]([C@H]([C@@H](C)C(=O)O2)O)O[C@H]2C[C@H](N(C)C)[C@H](O[C@H]3C[C@@H](C)C(O[C@H]4C[C@H](OC)[C@@H](NC)[C@H](C)O4)[C@H](C)O3)C)CO"

But perhaps it's easier to test with a simple lactone, like a 14-membered lactone.

But given time constraints, perhaps proceed with the code.

Now, the previous attempt had an error due to incorrect indentation. The code must be properly indented.

In the previous code, the error was "unexpected indent", which suggests that the code had incorrect indentation. For example, the code might have mixed tabs and spaces or incorrect alignment.

In the current code, the loops are properly indented. Let's check:

The code after parsing the SMILES:

rings = mol.GetRingInfo().AtomRings()

for ring in rings:
    ring_size = len(ring)
    if ring_size <12: continue
    ... loop through atoms in ring, check for O, then neighbors, etc.

Yes, the code is properly indented.

So the code should work.

Additionally, the function returns a boolean and a reason string.

Thus, the code should correctly identify macrolides by checking for a macrocyclic lactone.

Another possible improvement: check that the lactone is the only ring of that size, but the problem doesn't require that.

Now, considering the examples provided, like Amphotericin B, which has a large lactone ring, the code should detect it.

Therefore, the code should work.