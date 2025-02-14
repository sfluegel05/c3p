"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: CHEBI:35835 fatty alcohol

An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms. Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of a single alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) != 1:
        return False, f"Found {len(alcohol_matches)} alcohol groups, should be exactly 1"
    
    # Check allowed atom types
    allowed_atoms = set([6, 8, 1])  # C, O, H
    atom_nums = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    if not atom_nums.issubset(allowed_atoms):
        return False, "Molecule contains disallowed atom types"
    
    # Check number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3 or c_count > 27:
        return False, f"Found {c_count} carbon atoms, should be between 3 and 27"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 44 or mol_wt > 400:
        return False, "Molecular weight outside expected range for fatty alcohols"
    
    # Check saturation and branching
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable == c_count - 1:
        saturation = "saturated and unbranched"
    elif n_rotatable < c_count - 1:
        saturation = "unsaturated and/or branched"
    else:
        saturation = "highly branched and/or cyclic"
    
    return True, f"Aliphatic alcohol with {c_count} carbon atoms, {saturation}"