"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: CHEBI:33709 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino and carboxyl groups
    has_amino_group = any(atom.GetSymbol() == 'N' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 3 for atom in mol.GetAtoms())
    has_carboxyl_group = any(atom.GetSymbol() == 'C' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 3 and sum(1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'O') == 2 for atom in mol.GetAtoms())
    if not (has_amino_group and has_carboxyl_group):
        return False, "Missing amino and/or carboxyl group"

    # Check for a single chiral center with L configuration
    chiral_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and atom.GetHybridization() == Chem.HybridizationType.SP3]
    if len(chiral_centers) != 1:
        return False, "More than one chiral center or no chiral center found"

    chiral_center_idx = chiral_centers[0]
    chiral_center = mol.GetAtomWithIdx(chiral_center_idx)
    if chiral_center.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "Chiral center is not L configuration"

    # Check if the chiral center is attached to both amino and carboxyl groups
    amino_attached = False
    carboxyl_attached = False
    for neighbor in chiral_center.GetNeighbors():
        if neighbor.GetSymbol() == 'N' and sum(bond.GetBondTypeAsDouble() for bond in neighbor.GetBonds()) == 3:
            amino_attached = True
        elif neighbor.GetSymbol() == 'C' and sum(bond.GetBondTypeAsDouble() for bond in neighbor.GetBonds()) == 3 and sum(1 for nbr in neighbor.GetNeighbors() if nbr.GetSymbol() == 'O') == 2:
            carboxyl_attached = True
    if not (amino_attached and carboxyl_attached):
        return False, "Amino and carboxyl groups not attached to the chiral center"

    return True, "Molecule is a proteinogenic amino acid with L configuration"