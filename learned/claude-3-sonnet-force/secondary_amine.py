"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:35550 secondary amine
A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count nitrogen atoms
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(n_atoms) != 1:
        return False, "Must have exactly one nitrogen atom"
    
    # Get the nitrogen atom
    n_atom = n_atoms[0]
    
    # Check if nitrogen is aromatic
    if n_atom.GetIsAromatic():
        return False, "Nitrogen atom is aromatic, should be aliphatic"
    
    # Check for charges on nitrogen
    if n_atom.GetFormalCharge() != 0:
        return False, "Nitrogen atom has a formal charge"
    
    # Count substituents attached to nitrogen
    n_substituents = sum(1 for bond in n_atom.GetBonds() if bond.GetOtherAtom(n_atom).GetAtomicNum() not in [1, 6, 7, 8, 16])
    if n_substituents > 0:
        return False, "Nitrogen atom has additional substituents besides C, N, O, S"
    
    # Check if nitrogen has at least two substituents (C, N, O, S)
    n_subs = sum(1 for bond in n_atom.GetBonds() if bond.GetOtherAtom(n_atom).GetAtomicNum() in [6, 7, 8, 16])
    if n_subs < 2:
        return False, "Nitrogen atom has less than two substituents (C, N, O, S)"
    
    # Check if molecule is neutral
    mol_charge = rdMolDescriptors.CalcFormalCharge(mol)
    if mol_charge != 0:
        return False, f"Molecule has formal charge {mol_charge}, should be 0"
    
    return True, "Contains a secondary aliphatic amine group"