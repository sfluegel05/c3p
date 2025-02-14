"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.
    An N-acetyl amino acid is an amino acid where an acetyl group is attached to any nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acetyl amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for acetylated nitrogen
    # Matches any nitrogen atom attached to an acetyl group (N-C(=O)C)
    n_acetyl_pattern = Chem.MolFromSmarts("N-C(=O)C")
    acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    if not acetyl_matches:
        return False, "No N-acetyl group found"

    # SMARTS pattern for amino acid backbone (generalized)
    # Matches nitrogen connected to carbon connected to carboxyl group
    amino_acid_backbone_pattern = Chem.MolFromSmarts("N-C-C(=O)O")
    backbone_matches = mol.GetSubstructMatches(amino_acid_backbone_pattern)
    if not backbone_matches:
        return False, "No amino acid backbone found"
    if len(backbone_matches) != 1:
        return False, "More than one amino acid backbone found"

    # Check for peptide bonds (amide bonds between amino acids)
    # Peptide bond pattern: C(=O)-N-C-C(=O)O (excluding acetyl group)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N-C-C(=O)O")
    peptide_bond_matches = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomMapNum() == 0 and end_atom.GetAtomMapNum() == 0:
                sub_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()])
                if sub_mol.HasSubstructMatch(peptide_bond_pattern):
                    peptide_bond_matches.append(bond)
    if peptide_bond_matches:
        return False, "Peptide bond(s) found, molecule may be a peptide"

    # Ensure that the acetylated nitrogen is part of the molecule's nitrogen atoms
    n_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    acetyl_nitrogens = [match[0] for match in acetyl_matches]
    if not set(acetyl_nitrogens).intersection(n_indices):
        return False, "Acetyl group not attached to any nitrogen in the molecule"

    return True, "Molecule is an N-acetyl amino acid"