"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: Tertiary amine
A molecule is considered a tertiary amine if it has at least one neutral nitrogen atom
that is bonded to exactly three carbon atoms (and no hydrogens), representing substitution
of all three hydrogens of ammonia by hydrocarbyl groups.
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is defined as being formally derived from ammonia by replacing three
    hydrogen atoms by hydrocarbyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if a tertiary amine group is detected, False otherwise
        str: Explanation of the classification decision
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if atom is nitrogen (atomic number 7)
        if atom.GetAtomicNum() != 7:
            continue
        
        # We require a neutral nitrogen (formal charge 0)
        if atom.GetFormalCharge() != 0:
            continue
        
        # Get neighbors of the nitrogen atom
        neighbors = list(atom.GetNeighbors())
        
        # For a tertiary amine, there should be exactly 3 neighbors.
        if len(neighbors) != 3:
            continue
        
        # Further require that the nitrogen has no hydrogens (implicit or explicit)
        # This verifies that the nitrogen has been fully substituted.
        if atom.GetTotalNumHs() != 0:
            continue
        
        # Check that all three neighbors are carbon atoms.
        if all(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):
            return True, f"Found tertiary amine at atom index {atom.GetIdx()} with three carbon substituents."
    
    return False, "No tertiary amine group found in the molecule."

# Example usage (you can remove or comment out these lines in production):
if __name__ == "__main__":
    # Example: triethylamine: SMILES "CCN(CC)CC"
    test_smiles = "CCN(CC)CC"
    result, reason = is_tertiary_amine(test_smiles)
    print(result, reason)