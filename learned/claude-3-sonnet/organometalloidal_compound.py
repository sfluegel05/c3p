"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: CHEBI:50860 organometalloidal compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound has bonds between one or more metalloid atoms and 
    one or more carbon atoms of an organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define metalloids
    metalloids = {5: 'B', 14: 'Si', 32: 'Ge', 33: 'As', 51: 'Sb', 52: 'Te'}
    
    # Find metalloid atoms in molecule
    metalloid_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metalloids:
            metalloid_atoms.append(atom)
    
    if not metalloid_atoms:
        return False, "No metalloid atoms found"

    # Check for metalloid-carbon bonds
    for metalloid_atom in metalloid_atoms:
        metalloid_symbol = metalloids[metalloid_atom.GetAtomicNum()]
        
        # Get neighboring atoms
        neighbors = metalloid_atom.GetNeighbors()
        
        # Look for carbon neighbors
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Verify carbon is part of an organyl group by checking it has other carbon neighbors
                carbon_neighbors = neighbor.GetNeighbors()
                for carbon_neighbor in carbon_neighbors:
                    if carbon_neighbor.GetAtomicNum() == 6:  # Another carbon
                        return True, f"Found {metalloid_symbol}-C bond where C is part of an organyl group"
                
                # Special case: if the carbon has a double bond to oxygen (carbonyl)
                # or is part of a cyano group, it's still considered an organyl group
                for carbon_neighbor in carbon_neighbors:
                    if (carbon_neighbor.GetAtomicNum() == 8 and 
                        any(bond.GetBondType() == Chem.BondType.DOUBLE 
                            for bond in neighbor.GetBonds())):
                        return True, f"Found {metalloid_symbol}-C bond where C is part of a carbonyl group"
                    if carbon_neighbor.GetAtomicNum() == 7 and any(bond.GetBondType() == Chem.BondType.TRIPLE 
                                                                 for bond in neighbor.GetBonds()):
                        return True, f"Found {metalloid_symbol}-C bond where C is part of a cyano group"

    return False, "No metalloid-carbon bonds found where carbon is part of an organyl group"