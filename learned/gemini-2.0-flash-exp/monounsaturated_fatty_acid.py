"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    A monounsaturated fatty acid has one double or triple bond in the fatty acid chain,
    and a carboxylic acid group.
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a MUFA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Find the carbon atoms connected to the carboxyl group
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:
        return False, "No carboxylic acid group found"
    
    carboxyl_atom_idx = matches[0][0] #Get first atom of first match
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_atom_idx)

    # Get the carbon connected to carboxyl group (first neighbor not oxygen)
    first_chain_carbon = None
    for neighbor in carboxyl_carbon.GetNeighbors():
      if neighbor.GetAtomicNum() == 6:
        first_chain_carbon = neighbor
        break
    if first_chain_carbon is None:
        return False, "No carbon chain found connected to carboxyl group"


    def find_unsaturation(mol, current_atom_idx, visited_atoms, found_unsaturation):
        """Recursively search for a double or triple bond."""
        current_atom = mol.GetAtomWithIdx(current_atom_idx)
        if current_atom.GetAtomicNum() != 6:
          return False, False, visited_atoms
            
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in visited_atoms:
                continue
            
            bond = mol.GetBondBetweenAtoms(current_atom_idx, neighbor_idx)
            
            if bond and bond.GetBondType() == Chem.BondType.DOUBLE or bond.GetBondType() == Chem.BondType.TRIPLE:
                
                return True, True, visited_atoms

            
            if neighbor.GetAtomicNum() == 6:
                new_visited_atoms = visited_atoms.copy()
                new_visited_atoms.add(neighbor_idx)
                found_in_branch, found_overall, v_atoms = find_unsaturation(mol, neighbor_idx, new_visited_atoms, found_unsaturation)
                if found_in_branch:
                    return True, found_overall, v_atoms
        
        return False, found_unsaturation, visited_atoms

    # Initialize visited atoms with the carbon connected to COOH group
    visited_atoms = {first_chain_carbon.GetIdx()}

    has_unsaturation, found_overall, visited_atoms = find_unsaturation(mol, first_chain_carbon.GetIdx(), visited_atoms, False)

    if not has_unsaturation:
         return False, "No double or triple bond found in fatty acid chain"

    # Count the number of carbons in the chain
    c_count = len(visited_atoms)

    # Check for fatty acid chain length (at least 4 carbons in the chain, as smallest MUFA is butenoic acid)
    if c_count < 4:
        return False, f"Carbon chain length ({c_count}) is too short for a fatty acid"


    #Check if more than one double/triple bond
    num_double_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    num_triple_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C#C")))

    if (num_double_bonds + num_triple_bonds) != 1:
       return False, f"More than one double or triple bond found in molecule, should have exactly 1"

    return True, "Monounsaturated fatty acid identified"