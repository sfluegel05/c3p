"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:17164 primary amine
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is a compound where a nitrogen atom is bonded to one hydrocarbyl group and two hydrogens.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for primary amine: Nitrogen bonded to two hydrogens and one carbon
    primary_amine_pattern = Chem.MolFromSmarts("[#7H2][#6]")
    
    # Find matches in the molecule
    matches = mol.GetSubstructMatches(primary_amine_pattern)
    if not matches:
        return False, "No primary amine group found"
    
    # Check if nitrogen is only bonded to one carbon and two hydrogens
    for match in matches:
        nitrogen_idx = match[0]
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Check the total number of bonds (degree)
        if nitrogen.GetDegree() != 3:
            continue
        
        # Check that nitrogen is only bonded to one carbon and two hydrogens
        num_carbons = 0
        num_hydrogens = 0
        for neighbor in nitrogen.GetNeighbors():
            atomic_num = neighbor.GetAtomicNum()
            if atomic_num == 6:
                num_carbons += 1
            elif atomic_num == 1:
                num_hydrogens += 1
            else:
                break  # Nitrogen is bonded to another heteroatom
        else:
            if num_carbons == 1 and num_hydrogens == 2:
                return True, "Primary amine group found"
    
    return False, "Primary amine group not found or nitrogen bonded to other atoms"