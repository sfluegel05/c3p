"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile (any nitrile derived from an aliphatic compound)

An aliphatic nitrile is defined as any molecule that contains at least one nitrile
group (–C≡N) where the nitrile carbon is not aromatic, and its only substituent (the 
non‐N neighbor) is an aliphatic (carbon‐based) fragment that does not lead into any aromatic region.
"""

from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    
    The test requires that at least one nitrile group (C#N) in the molecule
    has a non‐nitrogen substituent that is itself a carbon (rejecting e.g. O–C#N as in sodium cyanate)
    and that this branch. when recursively examined, does not contain any aromatic atoms.
    (Also the substituent must have at least one hydrogen – a loose proxy for being an aliphatic fragment.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an aliphatic nitrile, else False.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the nitrile SMARTS pattern.
    nitrile_pattern = Chem.MolFromSmarts("[C;X2]#[N;X1]")
    if nitrile_pattern is None:
        return False, "Error creating nitrile pattern"
    
    # Get all nitrile matches.
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No nitrile group found in the molecule"
    
    # Recursive helper: check whether a branch (starting from an atom) contains any aromatic atoms.
    def branch_contains_aromatic(atom, banned_atom_idx, visited):
        """
        Recursively checks if from 'atom' (avoiding going back to banned_atom_idx)
        any aromatic atom is reached.
        """
        if atom.GetIdx() in visited:
            return False
        visited.add(atom.GetIdx())
        if atom.GetIsAromatic():
            return True
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == banned_atom_idx:
                continue
            if branch_contains_aromatic(nbr, banned_atom_idx, visited):
                return True
        return False

    # Loop over all nitrile substructure matches.
    for match in matches:
        # match[0] is the nitrile carbon; match[1] is the nitrile nitrogen.
        nitrile_c = mol.GetAtomWithIdx(match[0])
        nitrile_n = mol.GetAtomWithIdx(match[1])
        
        # Reject if the nitrile carbon is aromatic or not actually carbon.
        if nitrile_c.GetSymbol() != "C" or nitrile_c.GetIsAromatic():
            continue
        
        # Get the substituent(s) on the nitrile carbon that are not the nitrile nitrogen.
        substituents = [nbr for nbr in nitrile_c.GetNeighbors() if nbr.GetIdx() != nitrile_n.GetIdx()]
        if not substituents:
            continue  # no substituent found (unusual case)
        
        # We expect one substituent in a terminal nitrile.
        for sub in substituents:
            # The substituent must be a carbon (rejecting e.g. oxygen in [Na]OC#N)
            if sub.GetSymbol() != "C":
                continue
            # Optionally, require that the substituent has at least one hydrogen.
            if sub.GetTotalNumHs() < 1:
                continue
            
            # Recursively check that the branch (starting at the substituent) does not lead to any aromatic atom.
            # We ban going back to the nitrile carbon.
            if branch_contains_aromatic(sub, nitrile_c.GetIdx(), set()):
                continue
            
            # If we get here, we have a nitrile group with a substituent that appears
            # to be an aliphatic (carbon‐based) branch without any aromatic atoms.
            return True, "Contains a nitrile group attached to an exclusively aliphatic fragment"
    
    return False, "Nitrile group(s) found but none have an exclusively aliphatic substituent branch"


# Example usages (uncomment to test):
# print(is_aliphatic_nitrile("N#CCC#N"))  # malononitrile - expected True.
# print(is_aliphatic_nitrile("c1ccccc1C#N"))  # benzyl cyanide, but since the nitrile carbon is not directly attached to an aliphatic carbon, expected False.