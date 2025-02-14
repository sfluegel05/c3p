"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is a molecule with a single carboxylic acid group
    at the end of a carbon chain and can optionally have a single hydroxy group on the carbon chain
    (but not on the carboxyl carbon)

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise.
        str: Reason for the classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for allowed elements (C, H, O and D)
    allowed_elements = [1, 6, 8, 2]  # H, C, O, D
    for atom in mol.GetAtoms():
       if atom.GetAtomicNum() not in allowed_elements:
           return False, "Molecule contains disallowed elements"
    
    # Verify Carboxylic Acid group (must be exactly one with a protonated OH).
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
         return False, "Molecule does not contain exactly one carboxylic acid group"
    
    #label carboxyl carbon
    carbonyl_carbon = mol.GetAtomWithIdx(carboxylic_acid_matches[0][0])
    carbonyl_carbon.SetProp('carbonyl', 'True')

    # Check for saturation of carbons in the chain.
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Not a saturated molecule"
    
    # Check for branching by looking for carbons with more than two carbon neighbors, excluding the carboxyl carbon.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetSymbol() == 'C' and not atom.HasProp('carbonyl'): #carbon
           carbon_neighbors = 0
           for neighbor in atom.GetNeighbors():
               if neighbor.GetAtomicNum() == 6:
                    carbon_neighbors += 1
           if carbon_neighbors > 2:
               return False, "Branched chain detected"
    
    #check if the chain is cyclic, find the end of the chain (the carboxyl carbon) and check if a path exists back to this carbon.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.HasProp('carbonyl'):
          #get neighbors and check if we can get back to the original carbon without going through the carboxyl oxygen.
            neighbors = atom.GetNeighbors()
            if len(neighbors) < 2:
                 return False, "Invalid carboxyl group configuration"
            
            # check that there is exactly one neighbor which is a carbon atom
            carbon_neighbor_count = 0
            carbon_neighbor = None
            for n in neighbors:
                if n.GetAtomicNum() == 6:
                   carbon_neighbor_count += 1
                   carbon_neighbor = n
            
            if carbon_neighbor_count != 1:
               return False, "Invalid carboxyl group configuration"
            
            visited = {atom.GetIdx()}
            queue = [carbon_neighbor]
            
            while queue:
                current_atom = queue.pop(0)
                if current_atom.GetIdx() == atom.GetIdx():
                    return False, "Cyclic chain detected"
                
                visited.add(current_atom.GetIdx())
                
                for n in current_atom.GetNeighbors():
                    if n.GetAtomicNum() == 6 and n.GetIdx() not in visited:
                         queue.append(n)
                         
    # Check for hydroxy groups. Only allow 0 or 1 hydroxy group which is not bonded to the carboxyl carbon.
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    n_hydroxy = len(hydroxy_matches)

    if n_hydroxy > 1:
        return False, "Too many hydroxy groups"
    
    if n_hydroxy == 1:
       for match in hydroxy_matches:
           hydroxy_atom = mol.GetAtomWithIdx(match[0])
           #check is oxygen is connected to a carbon
           oxygen_neighbors = hydroxy_atom.GetNeighbors()
           if not any(neighbor.GetAtomicNum() == 6 for neighbor in oxygen_neighbors):
               return False, "Hydroxy group not connected to carbon chain"
           #check that the carbon attached to the hydroxy group is not the carbonyl carbon
           for neighbor in oxygen_neighbors:
               if neighbor.GetAtomicNum() == 6 and neighbor.HasProp('carbonyl') and neighbor.GetProp('carbonyl') == 'True':
                    return False, "Hydroxy group attached to carbonyl carbon"
                   
    
    # Check for other functional groups
    for atom in mol.GetAtoms():
       if atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 8:  # not hydrogen, carbon or oxygen. 
           return False, "Other functional groups detected"

    # Check number of carbon atoms (must be at least 4)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbon atoms"

    return True, "Straight-chain saturated fatty acid"