"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:50729 organoiodine compound
An organoiodine compound is a compound containing at least one carbon-iodine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for iodine atoms
    iodine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53]
    if not iodine_atoms:
        return False, "No iodine atoms found"

    # Check for carbon-iodine bonds
    c_i_bond_pattern = Chem.MolFromSmarts("[I;X1]-[C;X3]")
    c_i_bonds = mol.GetSubstructMatches(c_i_bond_pattern)
    if not c_i_bonds:
        return False, "No carbon-iodine bonds found"

    # Check that iodine is bonded directly to carbon in an organic moiety
    for bond in c_i_bonds:
        i_idx, c_idx = bond
        i_atom = mol.GetAtomWithIdx(i_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        
        # Check if carbon atom is part of an organic moiety
        if not any(nbr.GetAtomicNum() == 6 for nbr in c_atom.GetNeighbors()):
            return False, "Iodine not bonded to an organic carbon"
        
        # Check if iodine atom is bonded to only one carbon
        if len([bond for bond in mol.GetBondedAtoms(i_idx) if bond.GetAtomicNum() == 6]) > 1:
            return False, "Iodine bonded to multiple carbon atoms"

    return True, "Molecule contains at least one carbon-iodine bond in an organic moiety"