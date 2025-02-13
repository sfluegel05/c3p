"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:33838 amine

An amine is a compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Check for at least one N-C bond (nitrogen bound to carbon and not part of a ring)
    n_nc_bonds = sum(1 for bond in mol.GetBonds()
                     if bond.GetBeginAtom().GetAtomicNum() == 7 and bond.GetEndAtom().GetAtomicNum() == 6
                     and not bond.GetBeginAtom().IsInRing() and not bond.GetEndAtom().IsInRing())
    if n_nc_bonds == 0:
        return False, "No N-C bonds found (nitrogen not part of a ring)"

    return True, "Contains at least one nitrogen atom with N-C bond (not part of a ring)"