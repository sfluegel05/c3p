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
    
    # Check if molecule contains only C, H, N, O
    allowed_atoms = set([6, 1, 7, 8])  # C, H, N, O
    atom_nums = set([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains atoms other than C, H, N, O"
    
    # Check if molecule contains at least one nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro groups present"
    
    # Find the hydrocarbon skeleton
    hydrocarbon_frags = Chem.GetMolFrags(mol, aromHNAr=False, aromQAr=False, fragConnectedMol=False, sanitizeFrags=True)
    hydrocarbon_skeleton = None
    for frag in hydrocarbon_frags:
        if all(atom.GetAtomicNum() in [6, 1] for atom in frag.GetAtoms()):
            hydrocarbon_skeleton = frag
            break
    
    if hydrocarbon_skeleton is None:
        return False, "No hydrocarbon skeleton found"
    
    # Check if nitro groups are attached to the hydrocarbon skeleton
    nitro_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1]
    for nitro_atom in nitro_atoms:
        if not any(bond.GetBeginAtomIdx() == nitro_atom.GetIdx() for bond in hydrocarbon_skeleton.GetBonds()):
            return False, "Nitro group not attached to hydrocarbon skeleton"
    
    return True, "Molecule is a nitrohydrocarbon"