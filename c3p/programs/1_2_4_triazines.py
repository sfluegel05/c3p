"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: Chemical entities of the class 1,2,4-triazines
Definition: Any compound with a 1,2,4-triazine skeleton, in which nitrogen atoms replace carbon 
at positions 1, 2 and 4 of the core benzene ring structure.
Improvements over previous attempt:
  - Only consider six-membered rings containing exclusively C or N.
  - Ensure the ring contains exactly three Cs and three Ns.
  - Exclude fused rings by verifying that each atom in the candidate ring is only found in that ring.
  - Check all cyclic rotations for a N–N–C–N–C–C pattern.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine core based on its SMILES string.
    The core is defined as a six-membered ring made only of carbons and nitrogens (with exactly three
    of each) and that is not fused with another ring. Moreover, there is at least one cyclic permutation
    of the ring in which the atoms occur in the order (starting at some atom):
       position 0: nitrogen    (corresponds to benzene position 1)
       position 1: nitrogen    (benzene position 2)
       position 2: carbon      (benzene position 3)
       position 3: nitrogen    (benzene position 4)
       position 4: carbon      (benzene position 5)
       position 5: carbon      (benzene position 6)
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as having a 1,2,4-triazine skeleton, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule (includes sanitization)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information from the molecule
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # tuple of tuples of atom indices

    # Iterate over each ring
    for ring in all_rings:
        if len(ring) != 6:
            continue  # only interested in six-membered rings
        
        # Retrieve the atoms for the ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # First, ensure that every atom is either carbon or nitrogen.
        if any(atom.GetAtomicNum() not in (6, 7) for atom in ring_atoms):
            continue
        
        # Count nitrogens and carbons
        n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
        c_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        if n_count != 3 or c_count != 3:
            continue

        # Exclude rings that are fused with other rings.
        # For each atom in this ring, count how many rings (from all rings) it participates in.
        is_isolated = True
        for idx in ring:
            # Count how many rings (in all_rings) contain this atom index.
            count = sum(1 for r in all_rings if idx in r)
            if count > 1:
                is_isolated = False
                break
        if not is_isolated:
            continue
        
        # Try all cyclic rotations of the ring.
        # We need to see a permutation where indices (0, 1, 3) are nitrogen and (2, 4, 5) are carbon.
        for r in range(6):
            # Create rotated list of atoms
            rotated = [ring_atoms[(r + i) % 6] for i in range(6)]
            if (rotated[0].GetAtomicNum() == 7 and
                rotated[1].GetAtomicNum() == 7 and
                rotated[2].GetAtomicNum() == 6 and
                rotated[3].GetAtomicNum() == 7 and
                rotated[4].GetAtomicNum() == 6 and
                rotated[5].GetAtomicNum() == 6):
                return True, "Found 1,2,4-triazine ring pattern"
    
    return False, "No 1,2,4-triazine ring pattern found"