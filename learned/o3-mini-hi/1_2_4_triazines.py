"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: Chemical entities of the class 1,2,4-triazines
Definition: Any compound with a 1,2,4-triazine skeleton, in which nitrogen atoms replace carbon 
at positions 1, 2 and 4 of the core benzene ring structure.
"""

from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine core based on its SMILES string.
    The core is defined as a six-membered aromatic ring with exactly three nitrogen atoms and three carbon atoms,
    where the nitrogens occur in positions 1, 2 and 4 (when the benzene ring is numbered in sequence).
    
    For a six-membered ring with atoms numbered 0 to 5 in order, we require that there exists a rotation
    (cyclic permutation) such that:
    
      atom at index 0: nitrogen (position 1)
      atom at index 1: nitrogen (position 2)
      atom at index 2: carbon   (position 3)
      atom at index 3: nitrogen (position 4)
      atom at index 4: carbon   (position 5)
      atom at index 5: carbon   (position 6)
    
    This function examines all rings in the molecule looking for at least one match.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a 1,2,4-triazine skeleton, False otherwise.
        str : Reason for the classification.
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information from the molecule
    ring_info = mol.GetRingInfo()
    
    # Iterate over each ring in the molecule
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # only interested in six-membered rings
        
        # Retrieve atoms in the ring in the order given by ring indices
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # The ring must be aromatic to be considered a heteroaromatic core.
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        
        # Check that the ring contains exactly three nitrogens and three carbons
        n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
        c_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        
        if n_count != 3 or c_count != 3:
            continue
        
        # The next step is to validate the relative positions.
        # Because the order in the ring is arbitrary (cyclic), we try each rotation.
        for r in range(6):
            # Reorder the atoms so that our starting atom is at index 0.
            # Check positions (0, 1, and 3) should be nitrogen (atomic number 7),
            # and positions (2, 4, and 5) should be carbon (atomic number 6).
            pos0_is_N = ring_atoms[r].GetAtomicNum() == 7
            pos1_is_N = ring_atoms[(r+1) % 6].GetAtomicNum() == 7
            pos2_is_C = ring_atoms[(r+2) % 6].GetAtomicNum() == 6
            pos3_is_N = ring_atoms[(r+3) % 6].GetAtomicNum() == 7
            pos4_is_C = ring_atoms[(r+4) % 6].GetAtomicNum() == 6
            pos5_is_C = ring_atoms[(r+5) % 6].GetAtomicNum() == 6
            
            if pos0_is_N and pos1_is_N and pos2_is_C and pos3_is_N and pos4_is_C and pos5_is_C:
                return True, "Found 1,2,4-triazine ring pattern"
    
    return False, "No 1,2,4-triazine ring pattern found"