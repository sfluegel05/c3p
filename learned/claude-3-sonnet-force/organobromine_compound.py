"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:38241 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is a compound containing at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon-bromine bonds
    num_c_br_bonds = sum(bond.GetIsAromatic() and atom.GetAtomicNum() == 35 for bond in mol.GetBonds() for atom in [bond.GetBeginAtom(), bond.GetEndAtom()] if atom.GetAtomicNum() == 6)

    if num_c_br_bonds > 0:
        return True, f"Contains {num_c_br_bonds} carbon-bromine bond(s)"
    else:
        return False, "Does not contain any carbon-bromine bonds"