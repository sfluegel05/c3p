"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:33779 organic sulfide
Organic sulfides are compounds having the structure RSR (R =/= H)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide has the structure RSR (R != H).

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
    
    # Find sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    
    # Check if there is at least one sulfur atom
    if not sulfur_atoms:
        return False, "No sulfur atoms found"
    
    # Check if each sulfur atom is bonded to two non-hydrogen atoms
    for sulfur in sulfur_atoms:
        neighbors = [atom for bond in sulfur.GetBonds() for atom in [bond.GetBeginAtom(), bond.GetEndAtom()] if atom.GetAtomicNum() != 1]
        if len(neighbors) != 2:
            return False, f"Sulfur atom has {len(neighbors)} non-hydrogen neighbors, expected 2"
    
    # Check if all substituents on sulfur are not hydrogen
    for sulfur in sulfur_atoms:
        for bond in sulfur.GetBonds():
            neighbor = bond.GetOtherAtom(sulfur)
            if neighbor.GetAtomicNum() == 1:
                return False, "One of the substituents on sulfur is hydrogen"
    
    return True, "Contains the structure RSR where R is not hydrogen (organic sulfide)"