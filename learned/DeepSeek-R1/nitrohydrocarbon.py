"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
</think>

Okay, let's tackle this problem. The goal is to create a function called `is_nitrohydrocarbon` that checks if a given SMILES string represents a nitrohydrocarbon. From the definition, a nitrohydrocarbon is a hydrocarbon where one or more hydrogen atoms are replaced by nitro groups. Hydrocarbons are compounds made up of only carbon and hydrogen, so after substitution, the molecule should have nitro groups (NO2) attached but no other heteroatoms.

First, I need to parse the SMILES string into a molecule using RDKit. If the parsing fails, return False with an error message.

Next, check for the presence of nitro groups. A nitro group in SMILES is typically represented as [N+](=O)[O-]. Using RDKit's SMARTS pattern matching, I can search for this group. The SMARTS pattern for a nitro group would be "[N+](=O)[O-]". If there are no nitro groups, the molecule isn't a nitrohydrocarbon.

Then, ensure there are no other heteroatoms besides carbon, hydrogen, nitrogen in the nitro groups, and oxygen in the nitro groups. Wait, but hydrocarbons only have C and H. So after substitution with nitro groups, the allowed elements are C, H, N, and O (from the nitro). But wait, the nitro group is part of the substitution. So the original hydrocarbon only has C and H. The nitro group adds N and O. So the entire molecule should consist of C, H, N, and O, but the N and O must only be part of nitro groups.

Wait, but the nitro group is -NO2. So the molecule can have N and O only in the nitro groups. So any oxygen or nitrogen atoms must be part of a nitro group. Otherwise, if there's an oxygen in another functional group (like an alcohol, ether, etc.), it's not a nitrohydrocarbon.

So the steps would be:

1. Check that all nitrogen atoms are part of a nitro group.
2. Check that all oxygen atoms are part of a nitro group.
3. The rest of the atoms are carbon and hydrogen.

So, how to implement this?

For each atom in the molecule:
- If it's nitrogen, check that it's in a nitro group.
- If it's oxygen, check that it's in a nitro group.
- Carbon and hydrogen are allowed.

But how to verify that each N and O is part of a nitro group.

Alternative approach: After identifying all nitro groups, collect their atoms. Then, any N or O not in these groups would make the molecule invalid.

So steps:

- Find all nitro groups using SMARTS.
- Collect all N and O atoms in these groups.
- Iterate over all atoms in the molecule. For any atom that is N or O, check if it's in the collected set. If not, return False.

Also, the original hydrocarbon must have only C and H. So after substitution, the molecule can have C, H, N, O, but N and O must be exclusively in nitro groups.

Another thing: The nitro group is a substituent on the hydrocarbon. So the hydrocarbon structure should be such that the nitro groups are attached to carbon atoms. So the nitro groups should be connected to the carbon skeleton.

So the plan is:

1. Parse the SMILES into a molecule. Check validity.
2. Find all nitro groups ([N+](=O)[O-]).
3. If no nitro groups, return False.
4. Collect all N and O atoms in these nitro groups.
5. Check that all N and O atoms in the molecule are part of these nitro groups.
6. Ensure that the remaining atoms are only C and H.
7. Also, ensure that the nitro groups are attached to carbon atoms (since hydrocarbons have C and H, so substitution is on C).

Wait, but in the molecule, the nitro group is attached to a carbon. So each nitro group's nitrogen must be bonded to a carbon.

So perhaps, for each nitro group found, check that the nitrogen is bonded to a carbon.

Let me structure this:

First, find all nitro groups using SMARTS. The SMARTS pattern for a nitro group could be "[N+](=O)[O-]". However, in some cases, the nitro group might be represented without the charges, like "N(=O)=O". But in SMILES, the standard way is [N+](=O)[O-]. So using that pattern.

Then, for each match in the nitro groups, get the nitrogen atom. Check that it's bonded to a carbon.

But wait, the nitro group is attached to the hydrocarbon, which is a carbon chain. So the nitrogen in the nitro group must be connected to a carbon atom.

So, in the SMARTS pattern, perhaps we can include the connected carbon. For example, "C-[N+](=O)[O-]".

But maybe better to first find all nitro groups, then check their connectivity.

Alternatively, after finding all nitro groups, check that each nitro's nitrogen is bonded to a carbon.

So, code steps:

- Parse the molecule.
- Check for nitro groups.
- For each nitro group's nitrogen atom, check that it has a neighbor which is a carbon.
- Collect all N and O atoms in nitro groups.
- Iterate over all atoms in the molecule:
   - If atom is N or O, check if it's in the collected set. If not, invalid.
   - Otherwise, if it's not C or H, invalid (though hydrocarbons should only have C and H, but after substitution, nitro adds N and O).

Wait, but after substitution, the molecule has C, H, N, O. So any other elements (like S, Cl, etc.) would disqualify it.

So:

Check that all atoms are either C, H, N, or O.

Then, check that all N and O atoms are part of nitro groups.

So:

for atom in mol.GetAtoms():
    sym = atom.GetSymbol()
    if sym not in ['C', 'H', 'N', 'O']:
        return False, "Contains non-hydrocarbon elements besides nitro groups"
    if sym == 'N' or sym == 'O':
        if atom not in nitro_atoms:
            return False, "Nitrogen or oxygen not part of nitro group"

Now, how to collect nitro_atoms.

The SMARTS pattern for nitro group can be written as "[N+](=O)[O-]". Let's create a SMARTS pattern that matches the nitro group and captures the N and O atoms.

But when using GetSubstructMatches, each match is a tuple of atom indices that match the pattern. For the nitro group, the pattern would have three atoms: N, O (double bonded), and O- (single bonded). Wait, the nitro group is [N+](=O)[O-]. So the structure is N with a double bond to O and a single bond to O-.

So the SMARTS pattern would be "[N+:1](=[O:2])-[O-:3]". Then, each match would have three atoms: N, O (double bond), O- (single bond).

So when we find all matches, we can collect all N, O atoms from each match.

So:

nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
matches = mol.GetSubstructMatches(nitro_pattern)
nitro_atoms = set()
for match in matches:
    for atom_idx in match:
        nitro_atoms.add(mol.GetAtomWithIdx(atom_idx))

Then, any N or O not in nitro_atoms is invalid.

Putting this together.

Also, check that each nitro group is attached to a carbon. So for each nitro group's N atom, check that it's bonded to a carbon.

Wait, in the nitro group [N+](=O)[O-], the N is connected to the rest of the molecule. So in the SMILES string, the nitro group is attached to a carbon. For example, nitrobenzene is c1ccccc1[N+](=O)[O-]. So the N is bonded to a carbon in the benzene ring.

So for each nitro group's N atom, check that it has at least one neighbor which is a carbon.

So in the code:

for match in matches:
    # The nitro group's N is the first atom in the match (assuming the SMARTS is [N+](=O)[O-])
    n_idx = match[0]
    n_atom = mol.GetAtomWithIdx(n_idx)
    # Check neighbors of N atom
    has_carbon_neighbor = False
    for neighbor in n_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C':
            has_carbon_neighbor = True
            break
    if not has_carbon_neighbor:
        return False, "Nitro group not attached to carbon"

Putting all these steps together.

Now, testing the examples:

Take nitrobenzene: [O-][N+](=O)c1ccccc1. The code should find one nitro group. All N and O are part of nitro. The N is attached to a carbon (the benzene ring).

Another example: 2-nitropropane is CC(C)[N+]([O-])=O. Wait, the SMILES is written as CC(C)[N+](=O)[O-]. The nitro group is attached to the central carbon. So the code should detect that.

Now, possible edge cases:

- Molecules with other functional groups, like an alcohol (OH) would have an O not in a nitro group, so they should be rejected.
- Molecules with other nitrogen-containing groups (like amines) should be rejected.
- Deuterated compounds, like nitrobenzene-d5, which has deuterium instead of hydrogen. The code should still accept them as hydrocarbons (since D is an isotope of H, and in SMILES, H is represented as H or [H], and D as [2H] or [H]). RDKit treats isotopes as different atoms, but in terms of element, they are still H. Wait, no. The atomic number of deuterium is still 1 (same as H), but the isotope is different. So in RDKit, when checking the symbol, [2H] would have symbol 'H' with isotope 2. So in the code, when checking the elements, H includes all isotopes. So the code should still allow deuterated hydrocarbons.

Wait, in the code, when iterating over atoms, if the symbol is 'H' (including isotopes), then it's allowed. So deuterium would be considered as H in terms of element symbol.

But wait, in the code, the check is for atom.GetSymbol() in ['C', 'H', 'N', 'O']. For deuterium, the symbol is 'H' (since the element is hydrogen), even though it's an isotope. So the code would accept deuterated compounds as part of the hydrocarbon.

So nitrobenzene-d5 would pass, as all N and O are in nitro groups, and the rest are C and H (including deuterium).

Another example: a molecule with an oxygen in an ether group would have an O not in a nitro group, so it's rejected.

So the code should correctly identify nitrohydrocarbons and reject others.

Now, putting all this into code.

The function steps:

1. Parse the SMILES into a molecule. If invalid, return False.
2. Check that there is at least one nitro group.
3. Collect all atoms in nitro groups (N and O).
4. Check that all N and O atoms in the molecule are part of nitro groups.
5. Check that each nitro group's N is bonded to a carbon.
6. Ensure no other elements besides C, H, N, O.

Wait, but step 6 is redundant because step 4 checks that all N and O are in nitro groups, and the rest are C and H. So if any other elements (like S) are present, the atom.GetSymbol() would not be in ['C', 'H', 'N', 'O'], so step 4 would catch it.

Wait, no. The code first checks that all atoms are C, H, N, or O. Then, among those, N and O must be in nitro groups.

So in code:

After parsing the molecule, iterate over all atoms:

for atom in mol.GetAtoms():
    sym = atom.GetSymbol()
    if sym not in ['C', 'H', 'N', 'O']:
        return False, "Contains non-hydrocarbon elements besides nitro groups"

Then, check that all N and O are in nitro groups.

Also, check that each nitro's N is attached to C.

So the code:

from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon (a hydrocarbon with nitro groups replacing hydrogens).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for presence of nitro groups
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"
    
    # Collect all N and O atoms in nitro groups
    nitro_atoms = set()
    for match in nitro_matches:
        for atom_idx in match:
            nitro_atoms.add(mol.GetAtomWithIdx(atom_idx))
    
    # Check all atoms are C, H, N, O
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in ['C', 'H', 'N', 'O']:
            return False, f"Contains forbidden element: {sym}"
    
    # Check all N and O are in nitro groups
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in ['N', 'O'] and atom not in nitro_atoms:
            return False, f"{sym} atom not part of nitro group"
    
    # Check each nitro group's nitrogen is bonded to carbon
    for match in nitro_matches:
        n_idx = match[0]  # First atom in nitro pattern is the N
        n_atom = mol.GetAtomWithIdx(n_idx)
        has_carbon = False
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                has_carbon = True
                break
        if not has_carbon:
            return False, "Nitro group not attached to carbon"
    
    return True, "Hydrocarbon with nitro group substituents"

Testing this code with the examples:

- nitrobenzene: should return True. All N and O are in nitro, attached to C.
- 2-nitropropane: same.
- A molecule with an amine group: would have N not in nitro, so rejected.
- A molecule with an oxygen in an ester: O not in nitro, rejected.
- Deuterated compounds: H is allowed, so nitrobenzene-d5 passes.

This should cover the cases. Now, the previous attempt had an error due to incorrect indentation. Let's ensure the code is properly formatted.

Another thing: in the nitro_pattern, the SMARTS might need to be adjusted. Let me verify. The SMARTS [N+](=O)[O-] should match the nitro group. Let's test with nitrobenzene's SMILES: [O-][N+](=O)c1ccccc1. The nitro group is O-, N+ with double bond to O. The SMARTS should match this.

Yes, the nitro_pattern should correctly find the nitro group.

So this code should work.