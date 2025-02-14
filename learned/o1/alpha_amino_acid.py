"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is an amino acid in which the amino group is located on the carbon atom
    at the position alpha to the carboxy group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H1]')
    amino_group_pattern = Chem.MolFromSmarts('[NX3;H2,H1+]')
    
    # Find carboxylic acid groups
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    
    # Iterate over carboxylic acid groups
    for match in carboxylic_acid_matches:
        carboxy_carbon_idx = match[0]
        carboxy_carbon = mol.GetAtomWithIdx(carboxy_carbon_idx)
        
        # Find alpha carbon (carbon attached to carboxylic carbon that is not oxygen)
        neighbors = [atom for atom in carboxy_carbon.GetNeighbors() if atom.GetAtomicNum() != 8]
        if not neighbors:
            continue  # No alpha carbon found
        
        alpha_carbon = neighbors[0]
        
        # Check if alpha carbon is a carbon atom
        if alpha_carbon.GetAtomicNum() != 6:
            continue  # Alpha carbon is not carbon
        
        # Check for amino group attached to alpha carbon
        has_amino_group = False
        for alpha_neighbor in alpha_carbon.GetNeighbors():
            if alpha_neighbor.GetIdx() == carboxy_carbon_idx:
                continue  # Skip carboxylic carbon
            if alpha_neighbor.GetAtomicNum() == 7:
                # Check if nitrogen matches amino group pattern
                submol = Chem.PathToSubmol(mol, [alpha_neighbor.GetIdx()])
                if submol.HasSubstructMatch(amino_group_pattern):
                    has_amino_group = True
                    break
        if has_amino_group:
            return True, "Contains alpha-amino acid structure"
    
    return False, "No alpha-amino acid structure found"