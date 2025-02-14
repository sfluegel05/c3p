"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:35708 alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone contains a hydroxy group on the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carbonyl atoms
    carbonyl_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0 and atom.IsInRingSize(3)]

    for carbonyl_idx in carbonyl_atoms:
        # Get neighboring atoms of carbonyl
        neighbors = [atom.GetIdx() for atom in mol.GetAtomWithIdx(carbonyl_idx).GetNeighbors()]

        # Check if one neighbor is a carbon with hydrogen and the other neighbor is a carbon with oxygen
        alpha_carbon_idx = None
        hydroxy_carbon_idx = None
        for idx in neighbors:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 1:
                alpha_carbon_idx = idx
            elif atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 0:
                hydroxy_carbon_idx = idx

        # Check if we found both an alpha carbon and a hydroxy carbon
        if alpha_carbon_idx is not None and hydroxy_carbon_idx is not None:
            hydroxy_atom_idx = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(hydroxy_carbon_idx).GetNeighbors() if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1]
            if hydroxy_atom_idx:
                return True, "Contains a hydroxy group on the alpha-carbon relative to the C=O group"

    return False, "Does not contain an alpha-hydroxy ketone substructure"