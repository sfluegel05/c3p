"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid
Definition: Any dicarboxy acid carrying a hydroxy group on the carbon atom at position alpha to the carboxy group.
"""

from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is defined as any dicarboxylic acid where at least one
    of the carboxy groups has an adjacent (alpha) carbon bearing a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for carboxyl group: [CX3](=O)[OX2H]
    carboxyl_smarts = "[CX3](=O)[OX2H]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # Check that exactly two carboxyl groups are present
    if len(carboxyl_matches) != 2:
        return False, f"Molecule has {len(carboxyl_matches)} carboxylic acid groups (exactly 2 required)"
    
    # Add explicit hydrogens for reliable detection of hydroxy groups
    mol_h = Chem.AddHs(mol)
    
    # Initialize flag to determine if any carboxyl group possesses an alpha-OH group.
    alpha_hydroxy_found = False
    
    # Iterate over each carboxyl group match
    for match in carboxyl_matches:
        # According to our SMARTS match, the first atom is the carboxyl carbon.
        acid_carbon_idx = match[0]
        acid_carbon = mol_h.GetAtomWithIdx(acid_carbon_idx)
        
        # Identify the alpha carbon(s) attached to the carboxyl carbon.
        # We look for neighbors that are carbon and are not part of the carboxyl group itself.
        alpha_candidates = []
        for nbr in acid_carbon.GetNeighbors():
            # We expect the alpha carbon to be carbon
            if nbr.GetAtomicNum() == 6:
                alpha_candidates.append(nbr)
        
        # For each candidate alpha carbon, check if it has an -OH substituent.
        for alpha in alpha_candidates:
            # Check among the neighbors of alpha for an oxygen with at least one hydrogen
            for subnbr in alpha.GetNeighbors():
                if subnbr.GetAtomicNum() == 8:
                    # Check if this oxygen has any explicit hydrogen (using GetTotalNumHs())
                    if subnbr.GetTotalNumHs() >= 1:
                        alpha_hydroxy_found = True
                        break
            if alpha_hydroxy_found:
                break  # no need to check further if one alpha hydroxy group is identified
        if alpha_hydroxy_found:
            break
    
    if not alpha_hydroxy_found:
        return False, "No alpha carbon adjacent to a carboxyl group carries a hydroxyl group"
    
    return True, "Molecule is a 2-hydroxydicarboxylic acid (dicarboxylic acid with an alpha hydroxy substituent)"