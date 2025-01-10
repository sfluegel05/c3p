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
                # Further ensure neighbor depicts a lipid-like structure
                # Look for a chain of consecutive carbon atoms
                carbon_chain_length = 0
                visited = set()
                carbons_to_explore = [neighbor]

                while carbons_to_explore:
                    current_carbon = carbons_to_explore.pop()
                    if current_carbon.GetIdx() not in visited:
                        visited.add(current_carbon.GetIdx())
                        carbon_chain_length += 1
                        for further_neighbor in current_carbon.GetNeighbors():
                            if (further_neighbor.GetAtomicNum() == 6 and
                                further_neighbor.GetIdx() not in visited):
                                carbons_to_explore.append(further_neighbor)

                if carbon_chain_length > 8:  # Arbitrary threshold for a lipid-like chain
                    return True, "Contains sulfonic acid group bonded to a carbon with lipid-like structure"
        
    return False, "No carbon-sulfur bond found involving sulfonic acid"