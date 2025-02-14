"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: CHEBI:35741 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is a benzene ring with four chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count chlorine atoms
    n_chlorine = sum(atom.GetTotalNumHs(includeNeighbors=True) == 0 and atom.GetImplicitValence() == 1 for atom in mol.GetAtoms())
    
    # Tetrachlorobenzene must have exactly 4 chlorine atoms
    if n_chlorine != 4:
        return False, f"Found {n_chlorine} chlorine atoms, need exactly 4"
    
    # Check for benzene ring pattern
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"
    
    # Check if all chlorine atoms are attached to the benzene ring
    chlorine_pattern = Chem.MolFromSmarts("Cl")
    chlorine_atoms = mol.GetSubstructMatches(chlorine_pattern)
    
    benzene_atoms = mol.GetSubstructMatch(benzene_pattern)
    for cl_atom in chlorine_atoms:
        if not any(mol.GetBondBetweenAtoms(cl_atom, benzene_atom).GetIsAromatic() for benzene_atom in benzene_atoms):
            return False, "One or more chlorine atoms not attached to the benzene ring"
    
    return True, "Contains a benzene ring with 4 chlorine substituents"