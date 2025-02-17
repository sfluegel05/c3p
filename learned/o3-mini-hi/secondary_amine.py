"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine 
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
This function checks for any nitrogen atom in the molecule that has exactly two carbon (hydrocarbon) neighbors 
and one hydrogen (implicit or explicit), which conforms to a secondary amine.
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains a secondary amine functional group based on its SMILES string.
    A secondary amine is defined as an NH center that has two hydrocarbyl (carbon) groups attached.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains a secondary amine, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms in the molecule looking for a nitrogen (atomic number 7)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Get the total number of hydrogens (includes both implicit and explicit H atoms)
            hydrogen_count = atom.GetTotalNumHs()
            # Get the heavy atom neighbors (this does not count implicit hydrogens)
            heavy_neighbors = atom.GetNeighbors()
            # For a secondary amine the nitrogen should be attached to exactly 2 heavy atoms and 1 hydrogen
            if len(heavy_neighbors) == 2 and hydrogen_count == 1:
                # Check that both heavy substituents are carbon atoms (i.e. hydrocarbyl groups)
                if all(neighbor.GetAtomicNum() == 6 for neighbor in heavy_neighbors):
                    return True, "Found a nitrogen with exactly two carbon substituents and one hydrogen (secondary amine)"
    
    return False, "No secondary amine found: no nitrogen atom with exactly two carbon substituents and one hydrogen"