"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import numpy as np

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid has a hydroxyl group at the 16th carbon of a steroid core,
    where this hydroxy group has a beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a flexible steroid core SMARTS pattern
    steroid_core_smarts = "[C]1[C][C]2[C]3[C]([C]1)[C][C]4[C]3[C]([C]2)[C][C]5[C]4([C])[C]([C]5)"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)

    # Find substructure match for steroid core
    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if not core_match:
        return False, "No steroid core found"

    # Get the indices of core atoms
    core_atoms = list(core_match)

    # Identify C16: C16 is the carbon connected to the carbon at index 14 from the
    # SMARTS pattern, and the carbon at index 15 from the pattern
    c14_idx = core_atoms[13]
    c15_idx = core_atoms[14]

    # Find the neighbors of C14
    c14_atom = mol.GetAtomWithIdx(c14_idx)
    c14_neighbors = [neighbor.GetIdx() for neighbor in c14_atom.GetNeighbors()]
    
    # C16 is the common neighbor of C14 and C15
    c16_idx = None
    for neighbor in mol.GetAtomWithIdx(c15_idx).GetNeighbors():
        if neighbor.GetIdx() in c14_neighbors:
            c16_idx = neighbor.GetIdx()
            break

    if c16_idx is None:
          return False, "Cannot determine index of C16"

    # Get the C16 atom object
    c16_atom = mol.GetAtomWithIdx(c16_idx)
    
    #Check for a hydroxyl group attached to C16
    has_oh = False
    oh_oxygen_idx = None

    for neighbor in c16_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:  # Oxygen atom
          
            o_atom = neighbor
            has_oh = True
            oh_oxygen_idx = o_atom.GetIdx()
            break

    if not has_oh:
       return False, "No hydroxyl group found at C16"
        
    # Check for beta configuration
    # Method:
    # 1. get coordinates of C16, O, C13, and C17
    # 2. Construct a normal vector from C13 and C17, using the cross product
    # 3. Find the vector from C16 to O
    # 4. Project the C16-O vector onto the normal to find out the direction (z component)
    # If positive, is beta

    # Get coordinates and C13 and C17 indexes
    c13_idx = core_atoms[3]
    c17_idx = core_atoms[15]
    
    conf = mol.GetConformer()
    c13_coords = conf.GetAtomPosition(c13_idx)
    c17_coords = conf.GetAtomPosition(c17_idx)
    c16_coords = conf.GetAtomPosition(c16_idx)
    o_coords = conf.GetAtomPosition(oh_oxygen_idx)


    # Construct vectors
    v1 = np.array([c17_coords.x - c13_coords.x, c17_coords.y - c13_coords.y, c17_coords.z - c13_coords.z])
    
    # Get the index of C15 from the core_atoms and get C15's coordinates
    c15_idx = core_atoms[14]
    c15_coords = conf.GetAtomPosition(c15_idx)

    v2 = np.array([c15_coords.x - c13_coords.x, c15_coords.y - c13_coords.y, c15_coords.z - c13_coords.z])


    # Calculate the normal vector to the C13-C17 plane
    normal_vector = np.cross(v1, v2)
    
    # Calculate the C16-O vector
    c16_o_vector = np.array([o_coords.x - c16_coords.x, o_coords.y - c16_coords.y, o_coords.z - c16_coords.z])
    
    # Project onto the normal vector
    projection = np.dot(c16_o_vector, normal_vector)

    is_beta = projection > 0
    

    if not is_beta:
         return False, "Hydroxyl group is not beta-configured at C16"
    
    return True, "16beta-hydroxy steroid found"