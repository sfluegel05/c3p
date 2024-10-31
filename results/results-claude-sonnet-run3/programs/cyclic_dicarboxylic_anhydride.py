from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyclic_dicarboxylic_anhydride(smiles: str):
    """
    Determines if a molecule is a cyclic dicarboxylic anhydride.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cyclic dicarboxylic anhydride, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS pattern for cyclic anhydride: C(=O)OC(=O) in a ring
    patt = Chem.MolFromSmarts('[C;R]1(=O)[O;R][C;R](=O)[#6;R]~[#6;R]1')
    
    if not mol.HasSubstructMatch(patt):
        return False, "No cyclic anhydride pattern found"

    # Get matches and ring info
    matches = mol.GetSubstructMatches(patt)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    for match in matches:
        # Get the key atoms from the match
        c1, o1, o2, c2 = match[0:4]
        
        # Find the ring containing the anhydride group
        anhydride_ring = None
        for ring in rings:
            if c1 in ring and c2 in ring and o2 in ring:
                anhydride_ring = ring
                break
                
        if anhydride_ring is None:
            continue

        # Verify both carbonyls
        c1_atom = mol.GetAtomWithIdx(c1)
        c2_atom = mol.GetAtomWithIdx(c2)
        
        # Count oxygens connected to each carbonyl carbon
        c1_o_count = sum(1 for neighbor in c1_atom.GetNeighbors() if neighbor.GetAtomicNum() == 8)
        c2_o_count = sum(1 for neighbor in c2_atom.GetNeighbors() if neighbor.GetAtomicNum() == 8)

        # Count oxygens in the ring
        o_in_ring = sum(1 for atom_idx in anhydride_ring if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8)

        if (c1_o_count == 2 and c2_o_count == 2 and o_in_ring == 1):
            ring_size = len(anhydride_ring)
            return True, f"Contains cyclic anhydride in a {ring_size}-membered ring"

    return False, "No valid cyclic anhydride pattern found"
# Pr=1.0
# Recall=0.9090909090909091