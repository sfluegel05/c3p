"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: quaternary ammonium ion
A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded 
to nitrogen have been replaced with univalent (usually organyl) groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a quaternary ammonium ion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for nitrogen atoms
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not n_atoms:
        return False, "No nitrogen atoms found"
        
    # Check each nitrogen atom
    for n_atom in n_atoms:
        # Check formal charge
        if n_atom.GetFormalCharge() != 1:
            continue
            
        # Check number of bonds (should be 4 for quaternary ammonium)
        if len(n_atom.GetBonds()) != 4:
            continue
            
        # Check if all neighbors are carbons
        neighbors = n_atom.GetNeighbors()
        if not all(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):
            continue
            
        # Check if nitrogen has no hydrogens
        if n_atom.GetTotalNumHs() != 0:
            continue
            
        # If we get here, we've found a quaternary ammonium ion
        return True, "Contains nitrogen with +1 charge bonded to four carbon groups"
        
    return False, "No quaternary ammonium ion found"