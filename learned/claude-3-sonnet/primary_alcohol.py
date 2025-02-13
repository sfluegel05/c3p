"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: CHEBI:15734 primary alcohol

A primary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom which has either three hydrogen atoms attached to it or only one other carbon atom and two hydrogen atoms attached to it.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all hydroxy (-OH) groups
    hydroxy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1]
    
    # Check if at least one hydroxy group is present
    if not hydroxy_atoms:
        return False, "No hydroxy (-OH) group found"
    
    # Check each hydroxy group
    for hydroxy_atom in hydroxy_atoms:
        # Get the carbon atom attached to the hydroxy group
        carbon = hydroxy_atom.GetNeighbors()[0]
        
        # Check if the carbon is saturated (all bonds are single)
        if sum(bond.GetBondType() == Chem.BondType.SINGLE for bond in carbon.GetBonds()) != len(carbon.GetBonds()):
            continue
        
        # Check if the carbon has 3 hydrogens or 1 carbon + 2 hydrogens
        carbon_neighbors = carbon.GetNeighbors()
        hydrogen_count = sum(1 for neighbor in carbon_neighbors if neighbor.GetAtomicNum() == 1)
        carbon_count = sum(1 for neighbor in carbon_neighbors if neighbor.GetAtomicNum() == 6)
        
        if hydrogen_count == 3 or (hydrogen_count == 2 and carbon_count == 1):
            return True, "Contains a hydroxy (-OH) group attached to a saturated carbon atom with either 3 hydrogens or 1 carbon and 2 hydrogens"
    
    return False, "No primary alcohol group found"