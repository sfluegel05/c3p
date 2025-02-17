"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol (A compound that contains exactly two free/alcoholic –OH groups)
Definition:
    “Diol” means that the molecule has exactly two free (alcoholic) –OH groups attached to sp3 carbons.
    –OH groups that are part of carbonyl-based functionalities (e.g. carboxylic acids) are ignored
    because in those cases the oxygen is bound to an sp2 carbon.
    
The improved algorithm:
  1. Parse the SMILES string and add explicit hydrogens.
  2. Use a SMARTS pattern that encodes both the –OH group and its attachment to an sp3 carbon that is not
     double-bonded to oxygen. We use: "[C;sp3;!$([C]=[O])]-[O;X2H]".
  3. Count the unique –OH oxygen atoms that are found.
  4. Classify the molecule as a diol if and only if exactly two such free/alcoholic –OH groups are found.
     
Note:
  This approach still faces challenges on extreme and very complex structures. In those cases the “context”
  might be misinterpreted. Nevertheless it improves on the previous version by building the criteria into the match.
"""

from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is classified as a diol (contains exactly two free/alcoholic hydroxyl groups)
    based on its SMILES string.
    
    The function:
      - Converts the SMILES string into an RDKit molecule and adds explicit hydrogens.
      - Uses a SMARTS pattern that looks for an sp3 carbon (that is not part of a carbonyl) bearing an –OH.
      - Counts the unique –OH groups that satisfy the pattern.
      - Returns True if exactly two such groups are found.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diol (exactly two qualifying –OH groups), otherwise False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydroxyl hydrogens are visible
    mol = Chem.AddHs(mol)
    
    # SMARTS: an O atom with 2 connections and one hydrogen (alcoholic -OH), attached to a carbon atom that is:
    # a) sp3 and
    # b) does not have a double bond to oxygen (which would indicate a carbonyl)
    #
    # Note: The SMARTS "[C;sp3;!$([C]=[O])]-[O;X2H]" matches a substructure consisting of:
    # - a carbon atom (sp3, and not in a C=O) bonded (via a single bond) to
    # - an oxygen atom having exactly two connections (one to the carbon and one hydrogen).
    diol_smarts = "[C;sp3;!$([C]=[O])]-[O;X2H]"
    pattern = Chem.MolFromSmarts(diol_smarts)
    
    # Get the matches; each match is a tuple: (carbon_idx, oxygen_idx)
    matches = mol.GetSubstructMatches(pattern)
    
    # Collect unique oxygen indices from the matches.
    valid_oh_indices = set()
    for match in matches:
        # match[0] is the carbon, match[1] is the oxygen.
        oh_idx = match[1]
        valid_oh_indices.add(oh_idx)
    
    valid_oh_count = len(valid_oh_indices)
    
    if valid_oh_count == 2:
        return True, "Molecule contains exactly two free (alcoholic) hydroxyl groups attached to sp3 carbons and is classified as a diol."
    else:
        return False, f"Molecule contains {valid_oh_count} qualifying hydroxyl groups, which does not match the diol definition (exactly two required)."

# Example usage:
# test_smiles = "OCCCCCCCCCCCCO"  # Example: 1,12-dodecanediol
# result, reason = is_diol(test_smiles)
# print(result, reason)