"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:33838 amine

An amine is a compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least one nitrogen atom
    n_atoms = mol.GetNumAtoms()
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogen == 0:
        return False, "No nitrogen atoms found"

    # Check for at least one N-C bond (nitrogen bound to carbon)
    n_nc_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBeginAtom().GetAtomicNum() == 7 and bond.GetEndAtom().GetAtomicNum() == 6)
    if n_nc_bonds == 0:
        return False, "No N-C bonds found"

    # Check for at least one N-H bond (nitrogen bound to hydrogen)
    n_nh_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBeginAtom().GetAtomicNum() == 7 and bond.GetEndAtom().GetAtomicNum() == 1)
    if n_nh_bonds == 0:
        return False, "No N-H bonds found"

    # Check molecular weight - amines typically <300 Da
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt > 300:
        return False, "Molecular weight too high for amine"

    # Check for hydrogen deficiency - amines should be relatively saturated
    h_deficiency = mol.GetNumHeavyAtoms() - Descriptors.HeavyAtomCount(mol) + 1
    if h_deficiency > 4:
        return False, "Too many rings or double bonds for amine"

    return True, "Contains at least one nitrogen atom with N-C and N-H bonds"