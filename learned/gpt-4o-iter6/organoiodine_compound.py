"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine (C-I) bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of Iodine first for sanity check
    iodine_present = any(atom.GetAtomicNum() == 53 for atom in mol.GetAtoms())
    if not iodine_present:
        return False, "No iodine atoms in the molecule"

    # Check for presence of C-I bond using broad SMARTS pattern
    ci_bond_pattern = Chem.MolFromSmarts("[#6]-[#53]")  # Carbon to Iodine bond (single bond focused)
    if mol.HasSubstructMatch(ci_bond_pattern):
        return True, "Contains a carbon-iodine bond"

    # Check possibility of aromatic bonding or other pattern misses
    ci_aromatic_pattern = Chem.MolFromSmarts("[#6]:[#53]")  # Aromatic bond version
    if mol.HasSubstructMatch(ci_aromatic_pattern):
        return True, "Contains a carbon-iodine aromatic bond"

    return False, "No carbon-iodine bond found"