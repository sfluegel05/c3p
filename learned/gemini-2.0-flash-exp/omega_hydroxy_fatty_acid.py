"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (C(=O)O or C(O)=O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        carboxyl_pattern = Chem.MolFromSmarts("C(O)=O")
        if not mol.HasSubstructMatch(carboxyl_pattern):
            return False, "No carboxyl group found"

    # Check for at least one hydroxyl group (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
         return False, "No hydroxyl group found"
    
    # Find position of the carbonyl carbon in the COOH
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not carboxyl_matches:
        return False, "Could not determine carbonyl position"
    
    carbonyl_atom_idx = carboxyl_matches[0][0] 

    # Get atoms
    atoms = mol.GetAtoms()
    n_atoms = len(atoms)
    if n_atoms < 3:
        return False, "Too few atoms to be an omega-hydroxy fatty acid"
    
    # Count total number of carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Too few carbon atoms"
    

    #Find all terminal carbon atoms with no branching.
    terminal_carbon_pattern = Chem.MolFromSmarts("[CX4]([CX4])[#6]")
    terminal_carbon_matches = mol.GetSubstructMatches(terminal_carbon_pattern)

    valid_omega = False
    reason = "No omega-hydroxyl group found"
    
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    #Get the indices of the terminal carbon atoms
    terminal_carbons_indices = [match[2] for match in terminal_carbon_matches if match[2] != carbonyl_atom_idx ]


    for match in hydroxyl_matches:
            
            hydroxy_atom_idx = match[0]
            
            #Check if the hydroxy group is directly connected to one of the terminal carbon atoms
            hydroxy_atom = atoms[hydroxy_atom_idx]
            
            for neighbor in hydroxy_atom.GetNeighbors():
                if neighbor.GetIdx() in terminal_carbons_indices:
                  valid_omega = True
                  reason = "Omega hydroxyl group found"
                  break
            if valid_omega:
                break

    return valid_omega, reason