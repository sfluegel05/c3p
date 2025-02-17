"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic hydrocarbon moiety.
In our implementation, we require that:
    (1) The aldehyde carbon (CHO, with one hydrogen) is exocyclic (not part of a ring),
    (2) It is directly bonded to an aromatic carbon atom, and
    (3) At least one ring that contains that aromatic atom is composed entirely of carbon atoms.
Examples of valid compounds include piperonal, 2-hydroxy-1-naphthaldehyde, salicylaldehyde, etc.
Invalid (false‐positive) examples include heterocycle‐based aldehydes such as furan aldehydes.
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    
    An arenecarbaldehyde is defined as a molecule containing at least one aldehyde group (–CHO)
    where the carbonyl carbon is exocyclic and directly bonded to an aromatic carbon atom.
    Furthermore, the aromatic carbon atom must be part of at least one ring that is composed solely
    of carbon atoms (i.e. a true arene).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one valid arenecarbaldehyde group is detected, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS for an aldehyde group attached directly to an aromatic carbon.
    # [CX3H1](=O)[c] indicates:
    #    - a carbon atom with sp2 hybridization and exactly one hydrogen bound (i.e. CHO)
    #    - double-bonded to oxygen
    #    - and directly attached to an aromatic carbon ("[c]")
    aldehyde_smarts = "[CX3H1](=O)[c]"
    query = Chem.MolFromSmarts(aldehyde_smarts)
    if query is None:
        return False, "SMARTS query creation failed"
    
    matches = mol.GetSubstructMatches(query)
    if not matches:
        # Check if there is at least an aldehyde group that is not attached to an arene.
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3H1](=O)")):
            return False, "Aldehyde group present, but not exocyclicly attached to a pure aromatic (carbon-only) moiety."
        return False, "No aldehyde group detected."
    
    # Get ring information for the molecule (each ring is a tuple of atom indices).
    ringinfo = mol.GetRingInfo()
    rings = ringinfo.AtomRings()  # tuple of tuples
    
    # Check each match: we want one where:
    # a) The aldehyde carbon (first atom in match) is exocyclic.
    # b) The atom attached (second atom in match) is aromatic, is carbon (atomic number 6),
    #    and is in at least one ring that is composed solely of carbon atoms.
    for match in matches:
        aldehyde_carbon_idx = match[0]
        attached_atom_idx = match[1]
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_carbon_idx)
        attached_atom = mol.GetAtomWithIdx(attached_atom_idx)
        
        # Condition 1: Aldehyde carbon must be exocyclic
        if aldehyde_carbon.IsInRing():
            continue
        
        # Condition 2: The attached atom must be carbon (atomic number 6) and aromatic.
        if attached_atom.GetAtomicNum() != 6 or not attached_atom.GetIsAromatic():
            continue
        
        # Condition 3: Check that the attached aromatic carbon is in at least one ring 
        # where every aromatic atom is carbon.
        valid_ring_found = False
        for ring in rings:
            if attached_atom_idx in ring:
                # Check all atoms in the ring: if they are aromatic and not carbon, break.
                ring_valid = True
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetIsAromatic() and atom.GetAtomicNum() != 6:
                        ring_valid = False
                        break
                if ring_valid:
                    valid_ring_found = True
                    break
        if not valid_ring_found:
            continue
        
        # If we passed all three conditions, we have detected an arenecarbaldehyde group.
        return True, "Arenecarbaldehyde functional group detected (aldehyde directly attached exocyclically to an aromatic hydrocarbon moiety)."
    
    # If we found aldehyde groups but none with a successful attachment, return a negative result.
    return False, "Aldehyde group present, but not exocyclicly attached to a pure aromatic (carbon-only) moiety."


# Example usage (for testing purposes; uncomment to run tests):
# test_smiles = [
#     "[H]C(=O)c1ccc2OCOc2c1",         # piperonal (should be True)
#     "C12=CC=CC=C2C=CC(=C1C=O)O",       # 2-hydroxy-1-naphthaldehyde (should be True)
#     "[H]C(=O)c1c(OC)cc(\\C=C\\CCCC\\C=C\\C=C\\C)",  # false positive example
#     "[H]C(=O)C1=CC=C(C)O1",           # 5-methyl-2-furaldehyde (should be False)
#     "CC(=O)C"                        # aliphatic aldehyde (should be False)
# ]
# for s in test_smiles:
#     result, reason = is_arenecarbaldehyde(s)
#     print(s, "->", result, reason)