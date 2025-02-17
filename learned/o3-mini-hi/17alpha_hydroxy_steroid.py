"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17α-hydroxy steroid
Definition: The α-stereoisomer of 17-hydroxy steroid.
This code employs a heuristic:
  1. It parses the input SMILES.
  2. It checks that the molecule has at least 4 rings with ring sizes typical of a steroid nucleus (at least one 5-membered ring and at least three 6-membered rings).
  3. It looks for an –OH (hydroxyl) group.
  4. It then checks that at least one hydroxyl is attached to a chiral (stereochemically defined) carbon atom that is part of a ring.
Due to the complexity of assigning the exact C17 position and “α” orientation, this heuristic does not guarantee perfect classification.
If the molecule fails any test, a reason is returned.
If the molecule passes, the function returns True and a brief explanation.
If the task were too demanding, one might return (None, None).
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17α-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a 17α-hydroxy steroid, False otherwise 
        str: Explanation for the classification result
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings of the molecule.
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) < 4:
        return False, "Molecule does not have a fused four-ring system typical for steroids"
    
    # Count ring sizes. A typical steroid nucleus has three 6-membered rings and one 5-membered ring.
    ring_sizes = [len(r) for r in rings]
    if ring_sizes.count(6) < 3 or ring_sizes.count(5) < 1:
        return False, "Molecule does not contain the expected 3 six-membered and 1 five-membered rings typical of steroids"
    
    # Look for hydroxyl groups using a SMARTS for -OH.
    oh_smarts = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_smarts)
    if not oh_matches:
        return False, "No hydroxyl (–OH) group found"
    
    # Heuristic: Look for a hydroxyl attached to a ring carbon that is chiral.
    # (In a steroid, the 17-hydroxy is usually attached at a defined stereocenter.)
    candidate_found = False
    for match in oh_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Look at atoms bonded to the hydroxyl oxygen.
        for nb in o_atom.GetNeighbors():
            # Check that the neighbor is carbon and is in a ring.
            if nb.GetAtomicNum() == 6 and nb.IsInRing():
                # We require that the carbon has been assigned stereochemistry.
                # (Note: In many cases, the SMILES for steroids explicitly show stereochemistry.)
                if nb.HasProp('_ChiralityPossible') or nb.HasProp('_CIPCode'):
                    candidate_found = True
                    break
        if candidate_found:
            break
    if not candidate_found:
        return False, "No candidate hydroxyl-bearing chiral carbon found that could correspond to C17 with α-configuration"
    
    # If all tests pass, we assume the molecule is a 17α-hydroxy steroid.
    # (Again, note that a full classification would require alignment to a steroid template and explicit C17 assignment.)
    return True, "Molecule has a steroid-like four-ring system and a chiral hydroxyl group; heuristically classified as a 17α-hydroxy steroid"

# Example (uncomment to test):
# test_smiles = "CC(=O)[C@@]1(O)CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"  # Example steroid SMILES
# print(is_17alpha_hydroxy_steroid(test_smiles))