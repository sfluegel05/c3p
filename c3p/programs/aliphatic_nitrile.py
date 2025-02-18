"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile (any nitrile derived from an aliphatic compound)
"""

from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is defined as any nitrile where the -C≡N group is attached to an aliphatic (non‐aromatic) moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a nitrile group: a sp-hybridized carbon triple-bonded to a nitrogen.
    # [C;X2] ensures the carbon has two connections (one being the nitrile nitrogen and one additional bond).
    nitrile_pattern = Chem.MolFromSmarts("[C;X2]#[N;X1]")
    if nitrile_pattern is None:
        return False, "Error creating nitrile pattern"
    
    # Find all substructure matches for the nitrile group
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No nitrile group found in the molecule"
    
    # Check each nitrile group for aliphatic characteristics.
    # A nitrile is considered aliphatic if the nitrile carbon itself is not aromatic,
    # and the substituent (the neighbor other than the nitrile nitrogen) is also not aromatic.
    for match in matches:
        # match[0] is the nitrile carbon, match[1] is the nitrile nitrogen.
        nitrile_c = mol.GetAtomWithIdx(match[0])
        nitrile_n = mol.GetAtomWithIdx(match[1])
        
        # Check if the nitrile carbon is aromatic.
        if nitrile_c.GetIsAromatic():
            # Skip this nitrile group if the carbon is aromatic.
            continue
        
        # Identify the substituent attached to the nitrile carbon (ignore the nitrile nitrogen).
        # There can be only one neighbor aside from the nitrile nitrogen.
        substituent_found = False
        aliphatic_neighbor = True
        for neighbor in nitrile_c.GetNeighbors():
            if neighbor.GetIdx() == nitrile_n.GetIdx():
                continue  # Skip the nitrile nitrogen
            substituent_found = True
            # If this substituent is aromatic, then the nitrile is not derived from a purely aliphatic chain.
            if neighbor.GetIsAromatic():
                aliphatic_neighbor = False
                break
        
        # If there is no substituent (which is unusual for a nitrile) or if the substituent is aromatic,
        # then this nitrile group is not considered aliphatic.
        if not substituent_found:
            continue
        if not aliphatic_neighbor:
            continue
        
        # If we have at least one nitrile group that passes our tests, we consider the molecule an aliphatic nitrile.
        return True, "Contains a nitrile group attached to an aliphatic fragment"
    
    # If no nitrile groups met the aliphatic criteria:
    return False, "Nitrile group(s) found, but none are attached exclusively to aliphatic fragments"

# Example usages (you can uncomment these lines to test):
# print(is_aliphatic_nitrile("N#CCC#N"))  # malononitrile - should be True!
# print(is_aliphatic_nitrile("c1ccccc1C#N"))  # Benzyl cyanide (aromatic), should be False.