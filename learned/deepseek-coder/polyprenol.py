"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:26195 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a member of the class of prenols with the general formula H-[CH2C(Me)=CHCH2]nOH, 
    where the carbon skeleton is composed of more than one isoprene units.

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

    # Check for terminal hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[CH2][OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No terminal hydroxyl group found"

    # Define a more flexible isoprene unit pattern that accounts for different configurations
    isoprene_pattern = Chem.MolFromSmarts("[CH3][CH]=[CH][CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Need at least 2 isoprene units
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 2"

    # Check connectivity of isoprene units
    # Get all carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Find the terminal carbon with hydroxyl group
    terminal_carbon = None
    for atom in carbon_atoms:
        if atom.GetTotalNumHs() == 2 and len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]) == 1:
            terminal_carbon = atom
            break
    
    if terminal_carbon is None:
        return False, "Could not find terminal carbon with hydroxyl group"
    
    # Traverse the chain and count isoprene units
    visited = set()
    current_atom = terminal_carbon
    isoprene_count = 0
    prev_atom = None
    
    while current_atom is not None:
        visited.add(current_atom.GetIdx())
        # Check if current atom is part of an isoprene unit
        for match in isoprene_matches:
            if current_atom.GetIdx() in match:
                isoprene_count += 1
                break
        
        # Get next carbon in chain
        neighbors = [n for n in current_atom.GetNeighbors() 
                    if n.GetAtomicNum() == 6 and n.GetIdx() not in visited]
        
        # If multiple neighbors, choose the one that continues the chain
        if len(neighbors) > 1:
            # Prefer the neighbor that's part of an isoprene unit
            for n in neighbors:
                for match in isoprene_matches:
                    if n.GetIdx() in match:
                        current_atom = n
                        break
                else:
                    continue
                break
            else:
                current_atom = neighbors[0]
        elif len(neighbors) == 1:
            current_atom = neighbors[0]
        else:
            break
    
    # Need at least 2 isoprene units in the chain
    if isoprene_count < 2:
        return False, f"Found only {isoprene_count} isoprene units in the chain, need at least 2"

    # Check double bond pattern
    double_bonds = [bond for bond in mol.GetBonds() 
                   if bond.GetBondType() == Chem.BondType.DOUBLE]
    
    # Check that double bonds are in the expected positions
    # In polyprenols, double bonds should be between carbons 2 and 3 of each isoprene unit
    correct_double_bonds = 0
    for bond in double_bonds:
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        # Check if the double bond is between carbons 2 and 3 of an isoprene unit
        if (begin_atom.GetDegree() == 3 and end_atom.GetDegree() == 2) or \
           (begin_atom.GetDegree() == 2 and end_atom.GetDegree() == 3):
            correct_double_bonds += 1
    
    if correct_double_bonds < 2:
        return False, "Not enough correct double bonds for polyprenol structure"

    return True, "Contains terminal hydroxyl group and chain of more than one isoprene units with proper double bond pattern"