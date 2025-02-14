"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid has a hydroxyl group at carbon 17 with alpha stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid core with a 17-alpha-hydroxy group, using [CX4] for general substituted C
    # flexible ring system with [C], [CX3], [CX4] and 'c' for aromatic carbons.
    # C17 is defined as [C@] connected to O and must be in a ring and have 3 other connections
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]1[C][C][C]3[C]2[C][C][C]4[C]3[C@](O)([C,CX3,CX4])[C][C]4")

    if not mol.HasSubstructMatch(steroid_core_pattern):
          return False, "No steroid core with a 17-alpha-hydroxy group found."
    
    matches = mol.GetSubstructMatches(steroid_core_pattern)
    
    #If no matches for the substructure, return False
    if not matches:
        return False, "No steroid core match found"
    
    
    # we only take the first match, any match must be correct
    match = matches[0]

    # Get index of C17, the 16th atom in the SMARTS string
    c17_index_in_match = 15
    c17_atom_index = match[c17_index_in_match]

    c17_atom = mol.GetAtomWithIdx(c17_atom_index)
    
    
    # Get the oxygen directly connected to the c17, it should be an OH
    
    neighbors = c17_atom.GetNeighbors()
    
    found_hydroxy = False
    for neighbor in neighbors:
        if neighbor.GetAtomicNum() == 8:
            found_hydroxy = True
            break;
    
    if not found_hydroxy:
        return False, "No hydroxy group attached to C17"

    
    return True, "Contains a steroid core with a 17-alpha hydroxyl group"