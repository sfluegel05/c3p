"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is defined by an amide group derived from a fatty acid with significant carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of amide bond (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"
    
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Check for significant carbon chains linked to the amide bond
    significant_chain_length = 8  # Define threshold for chain as more than 8 carbon atoms

    for match in amide_matches:
        carbonyl_c_idx = match[0]  # Carbonyl Carbon
        nitrogen_idx = match[1]    # Nitrogen

        def count_carbon_chain(atom, exclude_idx=None):
            visited = set()
            carbon_chain = 0
            queue = [atom]
            while queue:
                current_atom = queue.pop(0)
                visited.add(current_atom.GetIdx())
                if current_atom.GetAtomicNum() == 6:  # Count carbon atoms only
                    carbon_chain += 1
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited and neighbor.GetIdx() != exclude_idx:
                        queue.append(neighbor)
            return carbon_chain

        # Count chain length from carbonyl carbon excluding nitrogen
        chain_from_carbonyl = count_carbon_chain(mol.GetAtomWithIdx(carbonyl_c_idx), nitrogen_idx)
        
        # If there's a significant carbon chain, classify as fatty amide
        if chain_from_carbonyl >= significant_chain_length:
            return True, "Contains amide bond derived from a significant carbon chain resembling fatty acid"
    
    return False, "Insufficient or no significant carbon chain beyond amide bond"