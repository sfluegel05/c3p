"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: CHEBI:35892 nitrile
A compound having the structure RC#N; thus a C-substituted derivative of hydrocyanic acid, HC#N.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains a cyano group (C#N) where the carbon is attached to another carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for C#N pattern
    cyano_pattern = Chem.MolFromSmarts("[C]#[N]")
    matches = mol.GetSubstructMatches(cyano_pattern)
    
    if not matches:
        return False, "No cyano group (C#N) found"
    
    # For each cyano group, verify it's attached to carbon (making it a nitrile)
    # and not to other atoms (which would make it a cyanide or other compound)
    for match in matches:
        c_atom = mol.GetAtomWithIdx(match[0])
        neighbors = c_atom.GetNeighbors()
        
        # Check if cyano carbon has exactly two neighbors (the N and one other atom)
        if len(neighbors) != 2:
            continue
            
        # Find the non-nitrogen neighbor
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() != 7:  # not nitrogen
                if neighbor.GetAtomicNum() == 6:  # is carbon
                    return True, "Contains C#N group with carbon attachment (nitrile group)"
                
    return False, "C#N group present but not properly connected to make a nitrile"