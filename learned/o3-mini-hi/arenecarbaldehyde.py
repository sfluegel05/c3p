"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
In our implementation, we require that:
    (1) The aldehyde carbon (–CHO) is exocyclic (i.e. not part of a ring),
    (2) It is directly attached to an aromatic carbon (atomic number 6 and flagged as aromatic),
    (3) That attached aromatic carbon is a member of at least one ring (of size 6 or more)
        that is composed solely of carbon atoms.
Examples include piperonal, salicylaldehyde, etc.
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    
    An arenecarbaldehyde contains an aldehyde group (–CHO) where:
      - The carbonyl carbon is exocyclic (not part of any ring).
      - It is connected to an aromatic carbon (atomic number 6, is aromatic).
      - That aromatic carbon is contained in at least one ring consisting exclusively of carbon atoms with length ≥ 6.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one valid arenecarbaldehyde group is detected, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS for a general aldehyde group (–CHO)
    aldehyde_query = Chem.MolFromSmarts("[CX3H1](=O)")
    if aldehyde_query is None:
        return False, "Failed to create aldehyde SMARTS"
    
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_query)
    if not aldehyde_matches:
        return False, "No aldehyde group detected."
    
    # Get ring info once from the molecule, to check ring membership later.
    rings = mol.GetRingInfo().AtomRings()
    
    # Loop over each aldehyde hit
    for match in aldehyde_matches:
        aldehyde_idx = match[0]  # This is the aldehyde carbon
        aldehyde_atom = mol.GetAtomWithIdx(aldehyde_idx)
        
        # Condition 1: Ensure the aldehyde carbon is exocyclic, i.e., not in any ring.
        if aldehyde_atom.IsInRing():
            # Skip this hit if the carbon is part of a ring.
            continue
        
        # For an aldehyde group, the carbon should be bonded to two atoms:
        # one (double-bonded) oxygen and one substituent.
        neighbors = aldehyde_atom.GetNeighbors()
        # Identify the oxygen double-bound partner and the substituent.
        substituent_atom = None
        for neighbor in neighbors:
            bond = mol.GetBondBetweenAtoms(aldehyde_idx, neighbor.GetIdx())
            # Identify the carbonyl oxygen based on bond type and atomic number
            if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                continue
            else:
                substituent_atom = neighbor
        if substituent_atom is None:
            continue  # no substituent found
        
        # Condition 2: Check that the substituent is an aromatic carbon.
        if substituent_atom.GetAtomicNum() != 6 or not substituent_atom.GetIsAromatic():
            continue
            
        # Condition 3: Check that the substituent carbon is part of at least one ring of size 6 or more
        # in which every atom is carbon.
        substituent_idx = substituent_atom.GetIdx()
        valid_ring_found = False
        for ring in rings:
            if substituent_idx in ring and len(ring) >= 6:
                # Ensure every atom in the ring is a carbon
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                    valid_ring_found = True
                    break
        if not valid_ring_found:
            continue
        
        # Passed all conditions: an arenecarbaldehyde group is identified.
        return True, "Arenecarbaldehyde functional group detected (aldehyde exocyclicly attached to an aromatic hydrocarbon moiety)."
    
    # If no matching hit passed our criteria, check if any aldehyde was found:
    if aldehyde_matches:
        return False, "Aldehyde group present, but not exocyclicly attached to an aromatic hydrocarbon moiety."
    return False, "No aldehyde group detected."