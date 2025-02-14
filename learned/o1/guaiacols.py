"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    
    A guaiacol is defined as any phenol carrying an additional methoxy substituent at the ortho-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for phenolic hydroxyl group
    # This pattern matches an aromatic carbon with an OH (phenol)
    phenol_pattern = Chem.MolFromSmarts('[cH]-[c]1[c][c][c][c][c]1[OH]')
    
    # Find all phenolic hydroxyl groups
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if not phenol_matches:
        return False, "No phenolic hydroxyl group found"
    
    # Define SMARTS pattern for methoxy group attached to an aromatic carbon
    methoxy_pattern = Chem.MolFromSmarts('[cH]-[c]([OCH3])[c]1[c][c][c][c][c]1')
    
    # Find all methoxy groups on aromatic carbons
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    if not methoxy_matches:
        return False, "No aromatic methoxy group found"
    
    # Check if any methoxy group is ortho to a phenolic OH
    for phenol_match in phenol_matches:
        phenol_carbon_idx = phenol_match[1]  # The aromatic carbon bonded to OH
        hydroxyl_oxygen_idx = phenol_match[-1]  # The OH oxygen atom
        
        # Get the neighboring carbons (ortho positions)
        phenol_carbon = mol.GetAtomWithIdx(phenol_carbon_idx)
        ortho_carbons = [nbr.GetIdx() for nbr in phenol_carbon.GetNeighbors() 
                         if nbr.GetIsAromatic() and nbr.GetAtomicNum() == 6]
        
        # Check if any ortho carbon has a methoxy group
        for ortho_idx in ortho_carbons:
            for methoxy_match in methoxy_matches:
                methoxy_carbon_idx = methoxy_match[1]  # Carbon atom bonded to OCH3
                if ortho_idx == methoxy_carbon_idx:
                    return True, "Molecule is a guaiacol (phenolic OH with ortho methoxy group)"
    
    return False, "No ortho methoxy group found adjacent to phenolic hydroxyl group"