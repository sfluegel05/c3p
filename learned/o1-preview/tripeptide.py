"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide that consists of three amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tripeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond pattern (C(=O)-N)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Find bonds corresponding to peptide bonds
    peptide_bond_bonds = []
    for match in matches:
        c_index = match[0]
        n_index = match[1]
        bond = mol.GetBondBetweenAtoms(c_index, n_index)
        if bond is not None:
            peptide_bond_bonds.append(bond.GetIdx())

    num_peptide_bonds = len(peptide_bond_bonds)

    # For tripeptide, there should be exactly 2 peptide bonds
    if num_peptide_bonds != 2:
        return False, f"Found {num_peptide_bonds} peptide bonds, expected 2 for tripeptide"

    # Break molecule at peptide bonds to get amino acid residues
    fragmented_mol = Chem.FragmentOnBonds(mol, peptide_bond_bonds, addDummies=False)
    fragments = Chem.GetMolFrags(fragmented_mol, asMols=True)

    num_fragments = len(fragments)
    if num_fragments != 3:
        return False, f"After breaking peptide bonds, found {num_fragments} fragments, expected 3"

    # Optionally, verify that each fragment resembles an amino acid (omitted for simplicity)

    return True, "Molecule is a tripeptide with three amino acid residues connected by peptide bonds"