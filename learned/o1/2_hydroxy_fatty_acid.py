"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:136568 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the alpha- or 2-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carboxylic acid group (protonated or deprotonated)
    carboxy_pattern = Chem.MolFromSmarts('[CX3](=O)[O;H1,-1]')
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # For each carboxylic acid carbon, check for 2-hydroxy group
    for match in carboxy_matches:
        c_idx = match[0]  # Carbonyl carbon atom index
        # Get the alpha carbon(s) connected to the carboxylic acid carbon
        carboxy_c_atom = mol.GetAtomWithIdx(c_idx)
        alpha_carbons = [nbr for nbr in carboxy_c_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]

        for alpha_c in alpha_carbons:
            # Check if alpha carbon has a hydroxy group attached
            has_hydroxy = False
            for nbr in alpha_c.GetNeighbors():
                if nbr.GetAtomicNum() == 8:  # oxygen atom
                    bond = mol.GetBondBetweenAtoms(alpha_c.GetIdx(), nbr.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        # Check if oxygen is a hydroxy group (attached hydrogen or implicit hydrogen)
                        if nbr.GetTotalDegree() == 1:
                            has_hydroxy = True
                            break
            if has_hydroxy:
                return True, "Contains 2-hydroxy group adjacent to carboxylic acid group"

    return False, "No 2-hydroxy group found adjacent to carboxylic acid group"