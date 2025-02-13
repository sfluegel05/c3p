"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a prenol with more than one isoprene units connected head-to-tail, ending with an -OH group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Normalize the molecule (add hydrogens)
    mol = Chem.AddHs(mol)

    # Define SMARTS pattern for isoprene unit connected head-to-tail
    isoprene_unit = Chem.MolFromSmarts("C(=C)C=C")  # Simplified pattern
    if isoprene_unit is None:
        return False, "Failed to create isoprene unit pattern"
    
    # Find all isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit)
    num_isoprene_units = len(isoprene_matches)
    
    if num_isoprene_units < 2:
        return False, f"Found {num_isoprene_units} isoprene units, need more than 1"
    
    # Check for hydroxyl group (-OH) at terminal position
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4;H2][OX2H]")  # Primary alcohol
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No terminal -OH group found"
    
    # Verify that isoprene units are connected head-to-tail
    # This is complex; we'll check that the isoprene units are connected in a chain
    # First, get atoms involved in isoprene units
    isoprene_atoms = set()
    for match in isoprene_matches:
        isoprene_atoms.update(match)
    
    # Check connectivity
    # Get the paths between the isoprene units
    is_connected = True
    isoprene_atoms_list = list(isoprene_atoms)
    for i in range(len(isoprene_atoms_list)-1):
        atom1 = mol.GetAtomWithIdx(isoprene_atoms_list[i])
        atom2 = mol.GetAtomWithIdx(isoprene_atoms_list[i+1])
        if not mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()):
            is_connected = False
            break
    
    if not is_connected:
        return False, "Isoprene units are not connected head-to-tail"
    
    return True, f"Contains {num_isoprene_units} isoprene units connected head-to-tail with terminal -OH group"