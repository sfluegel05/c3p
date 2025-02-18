"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: CHEBI:16669 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondStereo

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is an amino acid in which the amino group is located on the carbon atom
    at the position alpha to the carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find amino and carboxyl groups
    amino_groups = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() == 2]
    carboxyl_groups = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 and sum(bond.GetBondTypeAsDouble() == 2 for bond in atom.GetBonds()) == 2]

    if not amino_groups or not carboxyl_groups:
        return False, "No amino and/or carboxyl groups found"

    # Check if any amino group is alpha to a carboxyl group
    for amino_idx in amino_groups:
        amino_atom = mol.GetAtomWithIdx(amino_idx)
        for bond in amino_atom.GetBonds():
            alpha_atom = bond.GetOtherAtom(amino_atom)
            if alpha_atom.GetDegree() > 1:  # Avoid terminal atoms
                for alpha_bond in alpha_atom.GetBonds():
                    carboxyl_atom = alpha_bond.GetOtherAtom(alpha_atom)
                    if carboxyl_atom.GetIdx() in carboxyl_groups:
                        # Check stereochemistry
                        if alpha_bond.GetStereo() == BondStereo.STEREOANY or alpha_bond.GetStereo() == BondStereo.STEREONONE:
                            return True, "Contains an amino group at the alpha position to a carboxyl group"
                        elif alpha_bond.GetStereo() == BondStereo.STEREOZ:
                            if amino_atom.GetHybridization() == Chem.HybridizationType.SP3 and carboxyl_atom.GetHybridization() == Chem.HybridizationType.SP2:
                                return True, "Contains an amino group at the alpha position to a carboxyl group"
                        elif alpha_bond.GetStereo() == BondStereo.STEREOE:
                            if amino_atom.GetHybridization() == Chem.HybridizationType.SP2 and carboxyl_atom.GetHybridization() == Chem.HybridizationType.SP3:
                                return True, "Contains an amino group at the alpha position to a carboxyl group"

    return False, "No amino group found at the alpha position to a carboxyl group"