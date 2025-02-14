"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: CHEBI:35941 nitrohydrocarbon

A nitrohydrocarbon is a C-nitro compound that is a hydrocarbon in which one or more of the
hydrogens has been replaced by nitro groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.

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
    
    # Check if molecule contains any nitro groups
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro groups present"
    
    # Check if molecule is a hydrocarbon (contains only C, H, N, O)
    allowed_atoms = [6, 1, 7, 8]  # C, H, N, O
    atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if any(num not in allowed_atoms for num in atom_nums):
        return False, "Contains atoms other than C, H, N, O"
    
    # Check if nitro groups are attached to carbons
    nitro_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    for idx in nitro_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetTotalNumHs() > 0:
            return False, "Nitro group attached to hydrogen"
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                return False, "Nitro group not attached to carbon"
    
    return True, "Molecule is a nitrohydrocarbon"