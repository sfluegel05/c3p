"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
For example: piperonal ([H]C(=O)c1ccc2OCOc2c1) and salicylaldehyde ([H]C(=O)c1ccccc1O)
qualify as arenecarbaldehydes.

Our improved method uses a SMARTS query designed to capture exactly an aldehyde carbon
(that is, a carbon with one hydrogen and a double bond to oxygen) that is directly bonded
to an aromatic atom. Furthermore, we demand that the aldehyde carbon itself is not part of a ring,
so that the formyl group is exocyclic.
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as an aldehyde where the carbonyl carbon (CHO group)
    is exocyclic and is directly attached to an aromatic atom.

    We implement this by searching for the SMARTS pattern "[CX3H1](=O)[a]"
    (which ensures that the aldehyde carbon has exactly one bonded hydrogen and is
    double-bonded to oxygen, and is directly connected to an aromatic atom).
    We then additionally check that the aldehyde carbon is not itself part of a ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one valid arenecarbaldehyde group is detected, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS for an aldehyde group attached directly to an aromatic atom.
    # [CX3H1](=O)[a] means:
    #   - a tetra-coordinated (but sp2 hybridized) carbon with exactly one hydrogen (i.e. CHO)
    #   - double-bonded to an oxygen (=O)
    #   - and bonded to an aromatic atom ([a])
    aldehyde_smarts = "[CX3H1](=O)[a]"
    query = Chem.MolFromSmarts(aldehyde_smarts)
    if query is None:
        return False, "SMARTS query creation failed"
    
    matches = mol.GetSubstructMatches(query)
    
    # Check each match: we only want to validate those where the carbonyl carbon is
    # exocyclic (i.e. not part of a ring).
    for match in matches:
        # According to our query, the first atom in the match is the carbonyl carbon.
        aldehyde_carbon = mol.GetAtomWithIdx(match[0])
        if not aldehyde_carbon.IsInRing():
            return True, "Arenecarbaldehyde functional group detected (aldehyde directly attached to an aromatic moiety)."
    
    # If we found any aldehyde group at all, but none qualifies with an exocyclic attachment to an aromatic atom:
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3H1](=O)")):
        return False, "Aldehyde group present, but not exocyclicly attached to an aromatic moiety."
    
    return False, "No aldehyde group detected."

# Example usage (uncomment to test):
# print(is_arenecarbaldehyde("[H]C(=O)c1ccc2OCOc2c1"))  # Piperonal (expected True)
# print(is_arenecarbaldehyde("CC(=O)C"))                 # Aliphatic aldehyde (expected False)