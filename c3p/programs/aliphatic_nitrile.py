"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile (any nitrile derived from an aliphatic compound)

An aliphatic nitrile is defined as any molecule that contains at least one nitrile group (–C≡N)
where the nitrile carbon is not aromatic and where the nitrile carbon’s aliphatic substituent (i.e.
the neighbor other than the nitrile nitrogen) does not lead into any aromatic fragment.
"""

from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    A nitrile is considered aliphatic if the nitrile carbon is not aromatic,
    and its attached substituent (the non-nitrogen neighbor) does not lead to any aromatic atom.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an aliphatic nitrile, else False.
        str: Explanation for the classification
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a nitrile group.
    # [C;X2] ensures the carbon has exactly two bonds (one to the nitrile nitrogen and one more).
    nitrile_pattern = Chem.MolFromSmarts("[C;X2]#[N;X1]")
    if nitrile_pattern is None:
        return False, "Error creating nitrile pattern"
    
    # Find all nitrile substructure matches.
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No nitrile group found in the molecule"
    
    # Helper function to search the substituent branch for any aromatic atoms.
    def branch_contains_aromatic(atom, banned_atom_idx, visited):
        """
        Recursively check whether a branch from 'atom' contains any aromatic atoms.
        We will not traverse back to the nitrile carbon (banned_atom_idx) to avoid false hit.
        """
        if atom.GetIdx() in visited:
            return False  # already visited this atom, no aromatic found in this path yet.
        visited.add(atom.GetIdx())
        # If the current atom is aromatic, return True.
        if atom.GetIsAromatic():
            return True
        # Traverse neighbors (skipping the nitrile carbon which is banned).
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == banned_atom_idx:
                continue
            if branch_contains_aromatic(nbr, banned_atom_idx, visited):
                return True
        return False

    # For each nitrile group identified, perform tests.
    for match in matches:
        # match[0] is the nitrile carbon; match[1] is the nitrile nitrogen.
        nitrile_c = mol.GetAtomWithIdx(match[0])
        nitrile_n = mol.GetAtomWithIdx(match[1])
        
        # Reject if the nitrile carbon is aromatic.
        if nitrile_c.GetIsAromatic():
            continue
        
        # Identify the substituent on the nitrile carbon that is not the nitrile nitrogen.
        substituent_atoms = [nbr for nbr in nitrile_c.GetNeighbors() if nbr.GetIdx() != nitrile_n.GetIdx()]
        if not substituent_atoms:
            # Unusual nitrile if no substituent is found.
            continue
        # We typically expect only one substituent; if more than one, we check each.
        passed = False
        for sub in substituent_atoms:
            # First, if the immediate substituent is aromatic, skip it.
            if sub.GetIsAromatic():
                continue
            
            # Search the branch attached to the substituent (avoid going back to nitrile_c).
            if not branch_contains_aromatic(sub, nitrile_c.GetIdx(), set()):
                # We found a nitrile where the nitrile carbon and its non‐N neighbor (branch) are aliphatic.
                passed = True
                break
        
        if passed:
            return True, "Contains a nitrile group attached to an exclusively aliphatic fragment"
    
    return False, "Nitrile group(s) found, but none are attached exclusively to aliphatic fragments"

# Example usages:
# Uncomment to test the functionality:
# print(is_aliphatic_nitrile("N#CCC#N"))  # malononitrile - should be True.
# print(is_aliphatic_nitrile("c1ccccc1C#N"))  # Benzyl cyanide - should be False as the chain leads to aromatic ring.