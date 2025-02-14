"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:33121 organoiodine compound
An organoiodine compound is a compound containing at least one carbon-iodine bond.
"""
from rdkit import Chem

def has_carbon_iodine_bond(mol):
    """
    Checks if a molecule contains at least one carbon-iodine bond.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object

    Returns:
        bool: True if the molecule contains at least one carbon-iodine bond, False otherwise
    """
    for bond in mol.GetBonds():
        begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 53) or (begin_atom.GetAtomicNum() == 53 and end_atom.GetAtomicNum() == 6):
            return True
    return False

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine bond.

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
    
    # Check for carbon-iodine bonds
    if has_carbon_iodine_bond(mol):
        return True, "Contains at least one carbon-iodine bond"
    else:
        return False, "No carbon-iodine bonds found"