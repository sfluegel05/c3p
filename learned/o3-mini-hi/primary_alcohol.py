"""
Classifies: CHEBI:15734 primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Primary alcohol 
Definition: A primary alcohol is a compound in which a hydroxy group (-OH) is attached to a saturated carbon atom 
           that has either three hydrogen atoms attached (as in methanol) or only one other carbon atom and two hydrogen atoms attached (as in ethanol’s CH2OH group).
"""

from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule contains at least one primary alcohol group based on its SMILES string.
    
    A primary alcohol group is defined by a hydroxyl oxygen (–OH) attached to an sp3-hybridized carbon (the alpha carbon) 
    that carries either:
      - Three hydrogen atoms (making it CH3) OR
      - Two hydrogen atoms and one additional carbon neighbor (making it CH2OH).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if at least one primary alcohol group is detected, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms to find potential hydroxyl oxygens
    for atom in mol.GetAtoms():
        # Check if the atom is oxygen
        if atom.GetSymbol() != "O":
            continue

        # Check if oxygen has at least one hydrogen (implicit+explicit) attached. 
        # This indicates an -OH group rather than another oxygen-containing functionality.
        if atom.GetTotalNumHs() < 1:
            continue

        # Get neighbors of oxygen: we expect one of them to be the carbon bearing the -OH group.
        neighbors = atom.GetNeighbors()
        # Find neighbors that are carbon atoms
        carbon_neighbors = [nbr for nbr in neighbors if nbr.GetSymbol() == "C"]
        if len(carbon_neighbors) != 1:
            # If not exactly one carbon neighbor, this −OH is not in the desired alcohol environment.
            continue

        # Candidate is the carbon attached to this oxygen.
        cand = carbon_neighbors[0]
        # Check if the carbon is saturated (sp3 hybridized)
        if cand.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # Get the total number of hydrogens attached to the candidate carbon (implicit+explicit)
        h_count = cand.GetTotalNumHs()

        # Get heavy atom neighbors of the candidate carbon.
        # Note: In RDKit, GetNeighbors() returns only explicitly drawn heavy atoms.
        heavy_neighbors = cand.GetNeighbors()  
        # We know one of these heavy neighbors is the −OH oxygen.
        # The “other” neighbor (if present) should be a carbon for a CH2OH.
        other_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetIdx() != atom.GetIdx()]

        # Now check the count: either a CH3OH-type (methanol-like) or a CH2OH-type
        # In methanol (CH3OH): the candidate carbon has only the oxygen as heavy neighbor and 3 hydrogens.
        if h_count == 3 and len(heavy_neighbors) == 1:
            return True, "Found primary alcohol group (methanol-type, CH3OH)"
        # In other primary alcohols (CH2OH): candidate carbon has 2 hydrogens and exactly 2 heavy neighbors,
        # one of which is oxygen and the other must be carbon.
        elif h_count == 2 and len(heavy_neighbors) == 2:
            # Identify the non-oxygen neighbor
            nonO_neighbors = [nbr for nbr in other_neighbors if nbr.GetSymbol() == "C"]
            if len(nonO_neighbors) == 1:
                return True, "Found primary alcohol group (CH2OH-type)"
    
    return False, "No primary alcohol group found"
    
# (Optional) For test purposes you can include a __main__ block:
if __name__ == "__main__":
    # Example structures
    test_smiles = [
        "CC[C@@H](CO)[NH3+]",  # (2S)-2-aminobutan-1-ol: contains CH2OH
        "OCCCCCCCCCCCCCCCCCCCCCO",  # docosan-1-ol: primary alcohol at the end
        "OCc1ccco1",  # furfuryl alcohol: primary alcohol on a CH2OH attached to an aromatic moiety
        "NCCO",  # ethanolamine: has a CH2OH group (primary alcohol)
        "O=C(O)CO",  # glycolic acid: carboxylic acid with an -OH (not considered a simple alcohol here)
    ]
    
    for smi in test_smiles:
        result, reason = is_primary_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")