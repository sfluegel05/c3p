"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
Definition: An alkanesulfonate in which the carbon at position 1 is attached to R, 
which can represent hydrogens, a carbon chain, or other groups.
This program checks that there is a sulfonate group –S(=O)(=O)[O-] directly attached 
to an sp3 (non‐aromatic, acyclic) carbon.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    For our purposes we require that the molecule contains a –S(=O)(=O)[O-] group directly 
    attached to an sp3 (non-aromatic, acyclic) carbon. This carbon should represent the 
    terminal carbon (position 1) of an alkyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an alkanesulfonate oxoanion, False otherwise.
        str: Reason for the classification decision.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a sulfonate group linked to a non-aromatic carbon.
    # The pattern "[C;!a]S(=O)(=O)[O-]" matches:
    # • A non-aromatic carbon (C;!a)
    # • Directly attached to a sulfur atom with two =O groups and one [O-].
    sulfonate_pattern = Chem.MolFromSmarts("[C;!a]S(=O)(=O)[O-]")
    
    # Find all matches for the pattern.
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No sulfonate pattern (C–S(=O)(=O)[O-]) found"
    
    # Iterate over each match and apply extra filtering.
    # The match pattern returns indices: (carbon, sulfur, oxygen, oxygen, oxygen).
    for match in matches:
        # Get the carbon atom (first index in the match) and the sulfur (second index).
        carbon = mol.GetAtomWithIdx(match[0])
        sulfur = mol.GetAtomWithIdx(match[1])
        
        # Check 1: Carbon must be sp3.
        if carbon.GetHybridization() != rdchem.HybridizationType.SP3:
            continue  # Skip if the carbon is not sp3.
            
        # Check 2: The carbon should be acyclic (not part of a ring).
        if carbon.IsInRing():
            continue  # Skip if the carbon is in a ring.
        
        # Check 3: Aside from the bond to sulfur, none of the other heavy atom neighbors (atomic number > 1)
        # should be aromatic. This helps ensure the carbon is part of a simple alkyl chain.
        valid = True
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetIdx() == sulfur.GetIdx():
                continue  # Skip the sulfur bonded atom.
            # Consider only heavy atoms.
            if neighbor.GetAtomicNum() > 1 and neighbor.GetIsAromatic():
                valid = False
                break
        if not valid:
            continue  # Skip if an aromatic neighbor is found.
        
        # If we reached here, we have found a valid alkanesulfonate oxoanion moiety.
        return True, "Contains an alkanesulfonate oxoanion moiety attached to a terminal sp3 carbon"
    
    # If no matches pass all filters, then the sulfonate group is not in an appropriate alkane context.
    return False, "Found sulfonate group, but not attached to a terminal sp3 (acyclic, non-aromatic) carbon"
    
# (Optional) For testing purposes:
if __name__ == "__main__":
    test_smiles = [
        "OCCS([O-])(=O)=O",  # isethionate, should be True.
        "C1=CC=CS([O-])(=O)=O",  # An example with aromatic carbon attachment, should be False.
    ]
    for sm in test_smiles:
        result, reason = is_alkanesulfonate_oxoanion(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")