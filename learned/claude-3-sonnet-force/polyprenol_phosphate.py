"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: CHEBI:37412 polyprenol phosphate

A prenol phosphate resulting from the formal condensation of the terminal allylic 
hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group (P(O)(O)=O)
    phosphate_pattern = Chem.MolFromSmarts("P(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for linear isoprenoid chain (alternating double bonds and methyl groups)
    isoprenoid_pattern = Chem.MolFromSmarts("[C@H](C)=C[C@H](C)C=C")
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_pattern)

    if not isoprenoid_matches:
        return False, "No linear isoprenoid chain found"

    # Check if the isoprenoid chain is attached to the phosphate group
    for match in isoprenoid_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.HasSubstructMatch(phosphate_pattern):
                    return True, "Contains a linear isoprenoid chain attached to a phosphate group"

    return False, "Isoprenoid chain not attached to phosphate group"