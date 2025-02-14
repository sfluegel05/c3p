"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a molecule containing a benzene ring and four chlorine atoms.

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

    # Check for benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
         return False, "No benzene ring found"

    # Count all chlorine atoms in the molecule
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)

    if chlorine_count != 4:
        return False, f"Found {chlorine_count} chlorine atoms, need exactly 4"

    return True, "Molecule has a benzene ring and exactly four chlorine atoms"