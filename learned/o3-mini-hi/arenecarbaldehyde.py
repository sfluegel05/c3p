"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
In our implementation, we require that:
    (1) The aldehyde carbon (–CHO) is exocyclic (i.e. not part of a ring),
    (2) It is directly bonded to an aromatic carbon atom, and
    (3) That aromatic carbon is part of at least one ring (of size 6 or more) composed entirely of carbon atoms.
Examples of valid compounds include piperonal, 2-hydroxy-1-naphthaldehyde, salicylaldehyde, etc.
Invalid examples include aldehydes that are attached to heterocyclic moieties.
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    
    An arenecarbaldehyde is defined as a molecule containing at least one aldehyde group (–CHO)
    where:
      - The carbonyl carbon is exocyclic (not part of a ring).
      - It is directly attached to an aromatic carbon.
      - That aromatic carbon is a member of at least one ring (of size 6 or greater)
        that is composed exclusively of carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one valid arenecarbaldehyde group is detected, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS query for an aldehyde group attached directly to an aromatic carbon.
    # [CX3H1](=O)[c] means:
    #  - A carbon with three neighbors and one hydrogen (i.e. CHO)
    #  - Double-bonded to oxygen and then directly attached to an aromatic carbon.
    aldehyde_smarts = "[CX3H1](=O)[c]"
    query = Chem.MolFromSmarts(aldehyde_smarts)
    if query is None:
        return False, "SMARTS query creation failed"
    
    matches = mol.GetSubstructMatches(query)
    if not matches:
        # Check if there is an aldehyde group without the required attachment.
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3H1](=O)")):
            return False, "Aldehyde group present, but not exocyclicly attached to an aromatic hydrocarbon moiety."
        return False, "No aldehyde group detected."
    
    # Get ring information (each ring is represented by a tuple of atom indices).
    ringinfo = mol.GetRingInfo()
    rings = ringinfo.AtomRings()
    
    # Iterate over all matches to check if at least one meets our criteria.
    for match in matches:
        # According to our SMARTS, match[0] is the aldehyde carbon and match[1] is the aromatic carbon.
        aldehyde_carbon_idx = match[0]
        attached_atom_idx = match[1]
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_carbon_idx)
        attached_atom = mol.GetAtomWithIdx(attached_atom_idx)
        
        # Condition 1: The aldehyde carbon must be exocyclic.
        if aldehyde_carbon.IsInRing():
            continue
        
        # Condition 2: The attached atom must be an aromatic carbon.
        if attached_atom.GetAtomicNum() != 6 or not attached_atom.GetIsAromatic():
            continue
        
        # Condition 3: Check that the attached aromatic carbon participates in at least one ring
        # of size 6 or more that is composed solely of carbon atoms.
        found_valid_ring = False
        for ring in rings:
            if attached_atom_idx in ring and len(ring) >= 6:
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                    found_valid_ring = True
                    break
        if not found_valid_ring:
            continue
        
        # If all conditions are met, we have detected a valid arenecarbaldehyde group.
        return True, "Arenecarbaldehyde functional group detected (aldehyde is exocyclicly attached to an aromatic hydrocarbon moiety)."
    
    # If no match satisfies all criteria:
    return False, "Aldehyde group present, but not exocyclicly attached to an aromatic hydrocarbon moiety."
    
# Uncomment the following lines for testing with provided examples:
# test_smiles = [
#     "[H]C(=O)c1ccc2OCOc2c1",         # piperonal (expected True)
#     "C12=CC=CC=C2C=CC(=C1C=O)O",       # 2-hydroxy-1-naphthaldehyde (expected True)
#     "C1=C2CCN([C@](C2=C(C=3OCOC13)OC)([C@H](C=4C=CC(=C(C4C=O)OC)OC)OC(C)=O)[H])C", # 3-O-acetylpapaveroxine (True)
#     "O(C(=O)C(C)C)C1=C(OCC)C=C(C=C1)C([H])=O",                              # ethyl vanillin isobutyrate (True)
#     "[H]C(=O)c1ccc(O)cc1",           # salicylaldehyde (expected True)
#     "CC(=O)C",                      # aliphatic methyl ketone (expected False)
# ]
# for s in test_smiles:
#     result, reason = is_arenecarbaldehyde(s)
#     print(s, "->", result, reason)