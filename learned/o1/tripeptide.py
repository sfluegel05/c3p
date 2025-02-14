"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:16699 tripeptide
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide consisting of three amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amide bonds (peptide bonds)
    num_amide_bonds = rdMolDescriptors.CalcNumAmideBonds(mol)
    if num_amide_bonds != 2:
        return False, f"Expected 2 amide bonds for a tripeptide, found {num_amide_bonds}"

    # Identify peptide bonds using SMARTS pattern
    # Peptide bond: Carbonyl carbon double-bonded to oxygen and single-bonded to nitrogen
    peptide_bond_pattern = Chem.MolFromSmarts("[C;D2](=O)[NX3;H1,H2]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 2:
        return False, f"Expected 2 peptide bonds, found {len(peptide_bond_matches)}"

    # Fragment the molecule at peptide bonds
    bonds_to_break = []
    for match in peptide_bond_matches:
        carbon_idx = match[0]
        nitrogen_idx = match[1]
        bond = mol.GetBondBetweenAtoms(carbon_idx, nitrogen_idx)
        if bond is not None:
            bonds_to_break.append(bond.GetIdx())
    if len(bonds_to_break) != 2:
        return False, f"Expected to find 2 peptide bonds to break, found {len(bonds_to_break)}"

    # Create a new molecule with bonds broken
    fragmented_mol = Chem.FragmentOnBonds(mol, bonds_to_break, addDummies=False)

    # Get the individual fragments
    frags = Chem.GetMolFrags(fragmented_mol, asMols=True)
    if len(frags) != 3:
        return False, f"Expected 3 amino acid residues after fragmentation, found {len(frags)}"

    # Optionally, check if each fragment resembles an amino acid residue
    # Here, we can check if each fragment has an alpha carbon connected to an amine and a carboxyl group
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4H][CX3](=O)")
    for frag in frags:
        if not frag.HasSubstructMatch(amino_acid_pattern):
            return False, "One or more fragments do not resemble an amino acid residue"

    return True, "Molecule is a tripeptide with three amino acid residues connected by peptide bonds"