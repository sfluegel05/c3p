"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more hydrogens have been replaced by nitro groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitro group pattern -[N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro group found"
      
    # Check for carbon atom presence
    carbon_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "Does not contain a carbon backbone"

    # Evaluate presence of significant non-hydrocarbon elements but allow nitro derivs.
    allowed_elements = (1, 6, 7, 8, 16)  # H, C, N, O, S which are common in organonitros
    significant_non_hydrocarbon_atoms = [
        atom.GetAtomicNum() for atom in mol.GetAtoms() 
        if atom.GetAtomicNum() not in allowed_elements
    ]
    if len(significant_non_hydrocarbon_atoms) > 3:
        return False, "Contains too many non-hydrocarbon atoms unrelated to nitro groups"
        
    return True, "Molecule is a nitrohydrocarbon with one or more nitro groups attached to a carbon backbone"