"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for sulfonic acid group (S(=O)(=O)[O-])
    sulfonic_acid_pattern = Chem.MolFromSmarts("S(=O)(=O)O")
    matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
    
    if not matches:
        return False, "No sulfonic acid group found"
    
    # Check if there is a C-S bond involving the sulfur of the sulfonic acid group
    for match in matches:
        sulfur_idx = match[0]
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # Find carbon atoms bonded to this sulfur
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check for long carbon chain (basic lipid-like feature)
                # Add condition for lipid-like structure
                if len(list(neighbor.GetBonds())) >= 10:  # crude check for longer carbon chain
                    return True, "Contains sulfonic acid group bonded to a carbon with lipid-like structure"
        
    return False, "No carbon-sulfur bond found involving sulfonic acid"