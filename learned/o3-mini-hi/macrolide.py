"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
Definition: A macrocyclic lactone with a ring of twelve or more members derived from a polyketide.
A macrolide must contain an ester group (lactone) that is incorporated into a ring of at least 12 atoms.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is defined as a macrocyclic lactone with a ring of 12 or more members.
    
    The strategy is as follows:
    1. Parse the SMILES string into a molecule.
    2. Search for an ester substructure (C(=O)O).
    3. Iterate over rings (using the ring information from RDKit) and check if there is
       any ring with at least 12 atoms that contains both the ester carbon and the ester oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a macrolide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information for the molecule
    ring_info = mol.GetRingInfo().AtomRings()  # Returns a list of tuples, each tuple is a ring (atom indices)
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Define a SMARTS pattern for a simple ester group (lactone candidate)
    # The pattern "C(=O)O" will match three atoms: a carbonyl carbon, a doubly-bonded oxygen, and an ester oxygen.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found in the molecule"
    
    # Loop over all rings in the molecule
    for ring in ring_info:
        if len(ring) >= 12:
            # For each large ring, search if any ester match is contained (at least two key atoms of the ester)
            # We require that the carbonyl carbon and the ester oxygen (the first and final atoms in the match)
            # are both in the ring.
            for match in ester_matches:
                # match is a tuple of three atom indices corresponding to (C, =O, O)
                if match[0] in ring and match[2] in ring:
                    return True, f"Found macrocyclic lactone ring of size {len(ring)} containing an ester group"
    
    return False, "No macrocyclic lactone ring (ester embedded in ring with 12 or more atoms) found"
    
# Example usage (uncomment to test):
# test_smiles = "O1CCCCCCCCCC(OCCCCC1)=O"  # Example of a macrolactone ring (11 atoms may not qualify, adjust as needed)
# result, reason = is_macrolide(test_smiles)
# print(result, reason)