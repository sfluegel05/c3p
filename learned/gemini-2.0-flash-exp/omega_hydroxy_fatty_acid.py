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
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
         return False, "No hydroxyl group found"

    # Find position of the carbonyl carbon in the COOH
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    carbonyl_atom_idx = carboxyl_matches[0][0] if carboxyl_matches else None
    if carbonyl_atom_idx is None:
       return False, "Could not determine carbonyl position"


    # find all hydroxyl groups
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    
    # Get atoms
    atoms = mol.GetAtoms()
    n_atoms = len(atoms)
    if n_atoms < 3:
        return False, "Too few atoms to be an omega-hydroxy fatty acid"

    # Analyze the chain and find the terminal carbon atoms.
    # For each hydroxyl group, see if it is on the other end of a linear chain from the carbonyl group
    
    valid_omega = False
    reason = "No omega-hydroxyl group found"
    for match in hydroxyl_matches:
        hydroxyl_atom_idx = match[0]
        # Check each possible hydroxyl group to see if it is at the end of the chain opposite the carboxyl
        # The idea here is to try to find the number of bonds to the carboxyl
        
        
        #Check if the hydroxyl group is close to the carbonyl
        #Calculate the distance between carboxyl carbon and hydroxyl oxygen
        
        if hydroxyl_atom_idx == carbonyl_atom_idx:
           continue #same atom
        
        
        # Check if the hydroxyl group is on the end of chain opposite the carbonyl group
        
        
        #Get all paths from carbonyl to hydroxyl, if the shortest path passes through a 
        #large number of carbon atoms, then this is a valid omega hydroxy fatty acid.

        paths = Chem.GetShortestPaths(mol,carbonyl_atom_idx,hydroxyl_atom_idx)
        if not paths:
            continue #no path exists between carbonyl and hydroxyl
        
        
        #Check path to see if it is reasonable. 
        
        carbons_in_path = 0
        
        for path in paths:
            
            
           
            is_valid = True
            for i,atom_idx in enumerate(path):
               atom = atoms[atom_idx]
               if atom.GetAtomicNum() == 6:
                  carbons_in_path +=1
                  
               if i == 0 or i == len(path)-1:
                  if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 8 :
                     is_valid = False
                     break
               else:
                    if atom.GetAtomicNum() != 6:
                        is_valid = False
                        break
            if is_valid and carbons_in_path > 2:

                valid_omega = True
                reason = "Omega hydroxyl group found"
                break;

        if valid_omega:
            break;



    return valid_omega, reason