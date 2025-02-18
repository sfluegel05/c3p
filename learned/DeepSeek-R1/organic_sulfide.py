"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:48366 organic sulfide (RSR', R ≠ H)
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    Organic sulfides have the structure R-S-R' where R and R' are not hydrogen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate through all sulfur atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur
            # Check if sulfur has exactly two bonds and both are to non-hydrogen atoms
            if atom.GetDegree() == 2:
                neighbors = atom.GetNeighbors()
                # Ensure both neighbors are not hydrogens
                if all(n.GetAtomicNum() != 1 for n in neighbors):
                    return True, f"Contains sulfide group (R-S-R') with both R groups non-hydrogen"
    
    return False, "No sulfide group (R-S-R') with both R groups non-hydrogen found"