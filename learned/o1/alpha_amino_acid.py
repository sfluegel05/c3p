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
    
    # Define SMARTS patterns for carboxylic acid/carboxylate and amino groups
    carboxy_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]')
    amino_pattern = Chem.MolFromSmarts('[NX3;H2,H1+,H0]')
    
    # Find carboxylic acid/carboxylate groups
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid or carboxylate group found"
    
    # Iterate over carboxylic groups
    for carboxy_match in carboxy_matches:
        carboxy_carbon_idx = carboxy_match[0]
        carboxy_carbon = mol.GetAtomWithIdx(carboxy_carbon_idx)
        
        # Find alpha carbon (carbon attached to carboxylic carbon that is not oxygen)
        alpha_carbons = [neighbor for neighbor in carboxy_carbon.GetNeighbors() if neighbor.GetAtomicNum() == 6]
        if not alpha_carbons:
            continue  # No alpha carbon found
        
        for alpha_carbon in alpha_carbons:
            alpha_carbon_idx = alpha_carbon.GetIdx()
            
            # Check for amino group attached to alpha carbon
            amino_neighbors = [atom for atom in alpha_carbon.GetNeighbors() if atom.GetAtomicNum() == 7]
            if not amino_neighbors:
                continue  # No amino group attached to alpha carbon
            
            for amino_nitrogen in amino_neighbors:
                # Verify that the nitrogen matches the amino group pattern
                submol = Chem.PathToSubmol(mol, [amino_nitrogen.GetIdx()])
                if not submol.HasSubstructMatch(amino_pattern):
                    continue  # Not an amino group
                
                # Additional checks can be added here if necessary
                
                return True, "Contains alpha-amino acid structure"
    
    return False, "No alpha-amino acid structure found"