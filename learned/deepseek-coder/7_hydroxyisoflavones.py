"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: CHEBI:77064 7-hydroxyisoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone is an isoflavone with a hydroxy group at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible isoflavone core pattern with a hydroxy group at the 7-position
    # The pattern matches the benzopyran-4-one structure with a phenyl group at the 2-position
    # and a hydroxy group at the 7-position, allowing for substitutions on the phenyl ring.
    isoflavone_pattern = Chem.MolFromSmarts("[O]=[C]1[C]=[C][C]([OH])=[C][C]=[C]1-c2:c:c:c:c:c:2")
    
    # Check if the molecule matches the isoflavone core pattern with a hydroxy group at the 7-position
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Molecule does not match the 7-hydroxyisoflavone core structure"
    
    # Verify that the hydroxy group is at the 7-position
    # The 7-position is the carbon adjacent to the carbonyl group in the benzopyran-4-one core
    # We need to ensure that the hydroxy group is attached to this carbon
    # First, find the carbonyl oxygen
    carbonyl_oxygen = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:  # Oxygen with single bond (carbonyl)
            carbonyl_oxygen = atom
            break
    
    if carbonyl_oxygen is None:
        return False, "No carbonyl oxygen found"
    
    # Find the carbon adjacent to the carbonyl oxygen (7-position)
    seven_position_carbon = None
    for neighbor in carbonyl_oxygen.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon
            seven_position_carbon = neighbor
            break
    
    if seven_position_carbon is None:
        return False, "No carbon adjacent to carbonyl oxygen found"
    
    # Check if the 7-position carbon has a hydroxy group attached
    has_hydroxy = False
    for neighbor in seven_position_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:  # Hydroxy group
            has_hydroxy = True
            break
    
    if not has_hydroxy:
        return False, "No hydroxy group found at the 7-position"
    
    return True, "Molecule contains the 7-hydroxyisoflavone core structure with a hydroxy group at the 7-position"