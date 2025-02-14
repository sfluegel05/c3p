"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: CHEBI:15781 primary alcohol

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
    
    # Parse SMILES and handle tautomers
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return False, "Invalid SMILES string"
    tautomer = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
    
    # Find atoms with -OH group
    oh_atoms = [atom for atom in tautomer.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1]
    
    # Check if any of the -OH groups are attached to a primary carbon
    for oh_atom in oh_atoms:
        carbon = oh_atom.GetNeighbors()[0]
        if carbon.GetSymbol() == 'C':
            # Check if carbon has three hydrogen neighbors or one carbon neighbor and two hydrogens
            # Ignore other substituents like halogens, etc.
            hydrogen_count = sum(1 for neighbor in carbon.GetNeighbors() if neighbor.GetSymbol() == 'H')
            carbon_count = sum(1 for neighbor in carbon.GetNeighbors() if neighbor.GetSymbol() == 'C')
            if (hydrogen_count == 3) or (carbon_count == 1 and hydrogen_count == 2):
                return True, "Contains a primary alcohol group (-OH attached to a saturated carbon atom with either three hydrogen atoms or one other carbon atom and two hydrogen atoms)"
    
    # If no primary alcohol groups found
    return False, "Does not contain a primary alcohol group"