"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule that contains two amino acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for peptide bond (more specific)
    peptide_bond_smarts = '[C;D2](=O)[N;D2][C;H]'
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond is None:
        return False, "Invalid peptide bond SMARTS pattern"

    # Find all peptide bonds
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond)
    num_peptide_bonds = len(peptide_bond_matches)

    # A dipeptide should have exactly one peptide bond
    if num_peptide_bonds != 1:
        return False, f"Found {num_peptide_bonds} peptide bonds, need exactly 1"

    # Get the bond index of the peptide bond
    peptide_bond_atoms = peptide_bond_matches[0]
    peptide_bond_idx = mol.GetBondBetweenAtoms(peptide_bond_atoms[0], peptide_bond_atoms[1]).GetIdx()

    # Break the molecule at the peptide bond to get fragments
    frags = Chem.FragmentOnBonds(mol, [peptide_bond_idx], addDummies=True)
    frag_mols = Chem.GetMolFrags(frags, asMols=True)

    # Define the SMARTS pattern for an amino acid residue
    amino_acid_smarts = '[N;D2][C@@H]([C])[C](=O)'
    amino_acid = Chem.MolFromSmarts(amino_acid_smarts)
    if amino_acid is None:
        return False, "Invalid amino acid residue SMARTS pattern"

    # Count fragments that match an amino acid residue
    residue_count = 0
    for frag in frag_mols:
        if frag.HasSubstructMatch(amino_acid):
            residue_count += 1

    # Check if there are exactly two amino acid residues
    if residue_count == 2:
        return True, "Contains two amino acid residues connected by a peptide bond"
    else:
        return False, f"Found {residue_count} amino acid residues, need exactly 2"