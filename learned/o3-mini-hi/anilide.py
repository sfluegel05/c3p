"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide - Any aromatic amide obtained by acylation of aniline.
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined as an aromatic amide where the nitrogen of the amide bond
    is directly attached to an aromatic ring (i.e., derived from aniline).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an anilide, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a SMARTS pattern to find an amide functionality:
    # This pattern looks for a nitrogen (NX3) attached to a carbon (CX3) that is double-bonded to oxygen.
    # We expect the found pattern to represent part of an amide bond.
    amide_smarts = "[NX3][CX3](=O)"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    matches = mol.GetSubstructMatches(amide_query)
    
    if not matches:
        return False, "No amide functionality found"
    
    # For each amide match, check if the nitrogen is directly attached to an aromatic carbon.
    # The match tuple returns (nitrogen_index, carbon_index). We ignore the carbonyl oxygen in this pattern.
    for match in matches:
        n_idx = match[0]  # nitrogen atom of the amide
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check neighbors of the amide nitrogen that are not the carbonyl carbon
        aromatic_neighbor = False
        for neighbor in n_atom.GetNeighbors():
            # Skip the carbonyl group connection:
            if neighbor.GetAtomicNum() == 6 and neighbor.HasProp('_TripodAtom'):  # safeguard if used elsewhere; not required here
                pass  # ignore, not a robust label so we use index check instead
            # Ensure we are not considering the carbonyl carbon (which has a double bond to oxygen)
            if neighbor.GetIdx() == match[1]:
                continue
            # If the neighbor is carbon and is aromatic, then it likely represents the aniline part.
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                aromatic_neighbor = True
                break
        
        # If we find such an aromatic substitution on the amide nitrogen then we classify as anilide
        if aromatic_neighbor:
            return True, "Found amide group with the nitrogen directly attached to an aromatic ring (anilide)"
    
    return False, "Amide group(s) found but none exhibit a nitrogen attached directly to an aromatic ring"
    
# Example usage (you can test with one of the provided SMILES strings):
if __name__ == "__main__":
    test_smiles = "CC(=O)NC1=CC=C(C=C1)C(=S)NCC2=CC=CO2"  # one of the anilide examples
    result, reason = is_anilide(test_smiles)
    print("Is anilide:", result)
    print("Reason:", reason)