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
    
    # Define the SMARTS pattern for a sulfonic acid group with appropriate carbon linkage
    sulfonic_acid_smart = Chem.MolFromSmarts("S(=O)(=O)[O,N][CX4]")
    sulfonic_acid_matches = mol.GetSubstructMatches(sulfonic_acid_smart)
    
    if not sulfonic_acid_matches:
        return False, "No sulfonic acid group with a C linkage found"
    
    for match in sulfonic_acid_matches:
        sulfur_idx = match[0]
        carbon_link_idx = match[-1]
        
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        carbon_link_atom = mol.GetAtomWithIdx(carbon_link_idx)

        # Verify the carbon in question is part of a longer lipid-like chain
        visited = set()
        carbons_to_explore = [carbon_link_atom]
        carbon_chain_length = 0

        # Traverse the carbon chain
        while carbons_to_explore:
            current_carbon = carbons_to_explore.pop()
            if current_carbon.GetIdx() not in visited:
                visited.add(current_carbon.GetIdx())
                carbon_chain_length += 1
                for further_neighbor in current_carbon.GetNeighbors():
                    if further_neighbor.GetAtomicNum() == 6 and further_neighbor.GetIdx() not in visited and further_neighbor.GetDegree() > 1:
                        carbons_to_explore.append(further_neighbor)
        
        # Consider chains of minimum length representing a lipid
        if carbon_chain_length >= 10:
            return True, "Contains a sulfonic acid group bonded to a carbon that forms part of a lipid-like chain"
    
    return False, "No appropriate sulfonic acid group found forming a lipid-like structure"