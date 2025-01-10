"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: organohalogen compound
A compound containing at least one carbon-halogen bond (where X is a halogen atom)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define halogen atomic numbers
    halogens = {9: 'F', 17: 'Cl', 35: 'Br', 53: 'I'}
    
    # Look for halogen atoms in the molecule
    halogen_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in halogens:
            halogen_atoms.append(atom)
    
    if not halogen_atoms:
        return False, "No halogen atoms found"
    
    # Check if any halogen is bonded to carbon
    c_x_bonds = []  # List to store carbon-halogen bonds
    for halogen in halogen_atoms:
        # Get neighboring atoms
        neighbors = halogen.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atomic number
                halogen_symbol = halogens[halogen.GetAtomicNum()]
                c_x_bonds.append(f"C-{halogen_symbol}")
                
    if not c_x_bonds:
        return False, "No carbon-halogen bonds found"
    
    # Create detailed message about the types of C-X bonds found
    unique_bonds = set(c_x_bonds)
    bond_counts = {bond: c_x_bonds.count(bond) for bond in unique_bonds}
    bond_description = ", ".join(f"{count} {bond}" for bond, count in bond_counts.items())
    
    return True, f"Contains carbon-halogen bonds: {bond_description}"