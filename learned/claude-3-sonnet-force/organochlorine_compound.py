"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies chemical entities as organochlorine compounds (CHEBI:35486).
An organochlorine compound is a compound containing at least one carbon-chlorine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for direct carbon-chlorine bonds
    c_cl_bond_pattern = Chem.MolFromSmarts("[C]~[Cl]")
    has_c_cl_bond = mol.HasSubstructMatch(c_cl_bond_pattern)

    # Check for aromatic chlorine atoms
    aromatic_cl_pattern = Chem.MolFromSmarts("c-[Cl]")
    has_aromatic_cl = mol.HasSubstructMatch(aromatic_cl_pattern)

    # Check for chlorine atoms attached to heteroatoms
    heteroatom_cl_pattern = Chem.MolFromSmarts("[!#6;!#1]~[Cl]")
    has_heteroatom_cl = mol.HasSubstructMatch(heteroatom_cl_pattern)

    if has_c_cl_bond or has_aromatic_cl or has_heteroatom_cl:
        return True, "Contains at least one carbon-chlorine bond or chlorine atom attached to heteroatom or aromatic ring"
    else:
        return False, "No carbon-chlorine bonds or chlorine atoms attached to heteroatoms or aromatic rings found"