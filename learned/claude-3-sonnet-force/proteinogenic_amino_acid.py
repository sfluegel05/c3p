"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Generate tautomers and resonance structures
    tautomers = list(Chem.EnumerateTautomers(mol))
    all_mols = tautomers + [mol]

    for mol in all_mols:
        # Check for amino and carboxyl groups
        has_amino_group = any(atom.GetSymbol() == 'N' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 3 for atom in mol.GetAtoms())
        has_carboxyl_group = any(atom.GetSymbol() == 'C' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 3 and sum(1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'O') == 2 for atom in mol.GetAtoms())
        if not (has_amino_group and has_carboxyl_group):
            continue

        # Check for a single chiral center with L/S configuration
        chiral_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and atom.GetHybridization() == Chem.HybridizationType.SP3]
        if len(chiral_centers) != 1:
            continue

        chiral_center_idx = chiral_centers[0]
        chiral_center = mol.GetAtomWithIdx(chiral_center_idx)
        if chiral_center.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW and chiral_center.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            continue

        # Check if the chiral center is attached to both amino and carboxyl groups
        amino_attached = False
        carboxyl_attached = False
        for neighbor in chiral_center.GetNeighbors():
            if neighbor.GetSymbol() == 'N' and sum(bond.GetBondTypeAsDouble() for bond in neighbor.GetBonds()) == 3:
                amino_attached = True
            elif neighbor.GetSymbol() == 'C' and sum(bond.GetBondTypeAsDouble() for bond in neighbor.GetBonds()) == 3 and sum(1 for nbr in neighbor.GetNeighbors() if nbr.GetSymbol() == 'O') == 2:
                carboxyl_attached = True
        if amino_attached and carboxyl_attached:
            config = "L" if chiral_center.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW else "S"
            return True, f"Molecule is a proteinogenic amino acid with {config} configuration"

    return False, "Molecule is not a proteinogenic amino acid"