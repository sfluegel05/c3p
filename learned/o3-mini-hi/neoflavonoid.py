"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: Neoflavonoid 
Definition: Any 1-benzopyran with an aryl substituent at position 4.
Neoflavonoids are usually 2H-1-benzopyrans (i.e. chromenes) in which the non-fused carbon (position 4)
of the pyran ring bears an aryl (phenyl-like) substituent.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_neoflavonoid(smiles: str):
    """
    Determines if a given molecule (via its SMILES string) is a neoflavonoid.
    A neoflavonoid is defined as a 1-benzopyran (chromene) that has an aryl substituent (e.g., a phenyl ring)
    at the 4-position of the benzopyran core.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a neoflavonoid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS query for a 2H-1-benzopyran (chromene) with atom map numbers.
    # In a 2H-1-benzopyran the numbering is as follows:
    #   [O:1] -> position 1 (oxygen)
    #   [C:2] -> position 2 (fused with benzene)
    #   [C:3] -> position 3 (fused with benzene)
    #   [C:4] -> position 4 (non-fused, potential for aryl substituent)
    #   [c:5], [c:6], [c:7], [c:8] -> the benzene ring portion that fuses with the heterocycle.
    # The SMARTS below defines this bicyclic core.
    core_neoflav_smarts = "[O:1]1[C:2]=[C:3][C:4]1[c:5]1[c:6]c[c:7]c[c:8]1"
    core_neoflav_query = Chem.MolFromSmarts(core_neoflav_smarts)
    if core_neoflav_query is None:
        return False, "Error in SMARTS for benzopyran core"
    
    # Search for the benzopyran core in the molecule. There might be multiple matches.
    core_matches = mol.GetSubstructMatches(core_neoflav_query)
    if not core_matches:
        return False, "1-benzopyran (chromene) core not found"
    
    # Now check for an aryl substituent attached at the 4-position.
    # In our query the atom mapped as [C:4] is at index 3 in the match tuple.
    for match in core_matches:
        pos4_idx = match[3]
        pos4_atom = mol.GetAtomWithIdx(pos4_idx)
        
        # Look at neighbors of the position-4 atom that are not part of the benzopyran core.
        for neighbor in pos4_atom.GetNeighbors():
            if neighbor.GetIdx() in match:
                # neighbor is part of the core; skip it.
                continue
            # Check if the neighbor is aromatic.
            if not neighbor.GetIsAromatic():
                continue
            # As a heuristic, require that the neighbor is part of an aromatic ring (phenyl ring).
            # We check the molecule's ring info for a ring of size 6 that contains this neighbor and 
            # that all members are aromatic.
            rings = mol.GetRingInfo().AtomRings()
            for ring in rings:
                if neighbor.GetIdx() in ring and len(ring) == 6:
                    if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                        return True, "Found 1-benzopyran core with aryl substituent at position 4"
    
    return False, "No aryl substituent attached at position 4 of the 1-benzopyran core found"

# Example usage (you can remove these test lines when integrating into a larger package)
if __name__ == "__main__":
    # A sample SMILES string (one of the provided examples can be plugged in here)
    test_smiles = "CC1=CC(=CC2=C1C(CC3(O2)CC(NC(=S)N3)(C)C)C4=CC=C(C=C4)OC)O"
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)