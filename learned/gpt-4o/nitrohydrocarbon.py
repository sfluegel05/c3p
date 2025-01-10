"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more of the hydrogens has been replaced
    by nitro groups.

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
    
    # Define SMARTS pattern for nitro group attached to carbon
    nitro_pattern = Chem.MolFromSmarts("C[N+](=O)[O-]")
    
    # Search for nitro groups attached to carbon
    if mol.HasSubstructMatch(nitro_pattern):
        # Ensure it is a hydrocarbon: contains primarily C, H, and O
        valid_atoms = set([6, 1, 8, 7])  # Carbon, Hydrogen, Oxygen, Nitrogen
        unique_atoms = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
        
        if unique_atoms.issubset(valid_atoms):
            return True, "Contains nitro group(s) attached to hydrocarbon"

        return False, "Contains invalid atoms for hydrocarbon structure"

    return False, "No nitro groups attached to carbon found"