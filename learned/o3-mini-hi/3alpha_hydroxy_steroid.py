"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI: 3α-hydroxy steroid (heuristic)
A 3α-hydroxy steroid is defined as a steroid having a fused tetracyclic core (typically three six‐membered rings and one five‐membered ring),
with a hydroxyl (-OH) substituent on one of the six‐membered rings of that core.
This heuristic verifies that the molecule has multiple fused rings and then checks that at least one hydroxyl group is attached
to a carbon that is in a six-membered ring and that is part of at least two rings (indicative of fusion).
"""

from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3α-hydroxy steroid based on its SMILES string using an improved heuristic.
    
    The heuristic performs these steps:
      1. Parse the SMILES string.
      2. Check that the molecule has at least 4 rings (a typical steroid has 4 fused rings).
      3. Identify all hydroxyl groups (i.e. –OH groups: oxygen with a hydrogen).
      4. For each hydroxyl group, check if its attached carbon is in any six-membered ring 
         (by examining the ring membership of that atom) and is part of at least 2 rings (as a proxy for being fused).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3α-hydroxy steroid; False otherwise.
        str: A message explaining the reasoning behind the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    total_rings = ring_info.NumRings()
    # Typical steroids have a fused tetracyclic system (4 rings)
    if total_rings < 4:
        return False, f"Insufficient rings found: {total_rings} rings, while steroids typically have 4 fused rings"
    
    # Get a list of all ring atom index tuples
    rings = ring_info.AtomRings()
    
    # Create a lookup: for each atom index, list the sizes of rings in which it participates.
    atom_ring_sizes = {}
    for ring in rings:
        ring_size = len(ring)
        for idx in ring:
            atom_ring_sizes.setdefault(idx, []).append(ring_size)
    
    # Search for all hydroxyl groups (oxygen with one bond to carbon and one hydrogen)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    if oh_pattern is None:
        return False, "Error generating the hydroxyl pattern"
    
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if not oh_matches:
        return False, "No hydroxyl (-OH) groups found in the molecule"
    
    # For each hydroxyl group found, check if its neighbor is a carbon that is in a six-membered ring AND in at least 2 rings
    for match in oh_matches:
        o_idx = match[0]  # index of the -OH oxygen
        o_atom = mol.GetAtomWithIdx(o_idx)
        # inspect neighbors; typically the oxygen is connected to one carbon atom
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue
            c_idx = neighbor.GetIdx()
            # Check if this carbon is part of any six-membered ring
            # and that it participates in 2 or more rings (a simple proxy for being in a fused ring system)
            if c_idx in atom_ring_sizes:
                # Find if one of the rings is six-membered
                in_six = any(size == 6 for size in atom_ring_sizes[c_idx])
                if in_six and (len(atom_ring_sizes[c_idx]) >= 2):
                    return True, ("Molecule has at least 4 rings and contains a hydroxyl group attached to a six-membered ring carbon "
                                  "that is in a fused environment (>=2 rings), consistent with a 3α-hydroxy steroid")
    
    return False, "No hydroxyl group found on a six-membered ring carbon in a fused ring system of the molecule"


# Example usage:
if __name__ == "__main__":
    # Test on one provided example:
    test_smiles = "C1[C@@]2([C@@]([C@@]3([C@](C[C@H](O)CC3)(C1)[H])C)(CC[C@@]4([C@H](CC[C@@]24[H])OS(O)(=O)=O)C)[H])[H]"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)