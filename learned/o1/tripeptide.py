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

    # Define a more specific peptide bond pattern
    # Peptide bond: O=C-NH, where the carbonyl carbon is connected to an amide nitrogen
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)-N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    # Define alpha amino acid pattern
    # Alpha amino acid: N-C-C(=O), where N is attached to C(alpha), which is attached to a carbonyl group
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1][CX4H][CX3](=O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    num_amino_acids = len(amino_acid_matches)

    # For tripeptides, we expect 3 amino acid residues
    if num_amino_acids != 3:
        return False, f"Found {num_amino_acids} amino acid residues, need exactly 3 for tripeptide"

    # Ensure that the amino acids are connected via peptide bonds
    # Check the connectivity of the peptide bonds
    valid_connections = 0
    for match in peptide_bond_matches:
        # Get atoms involved in the peptide bond
        carbonyl_c = mol.GetAtomWithIdx(match[0])
        amide_n = mol.GetAtomWithIdx(match[1])

        # Check if the amide nitrogen is connected to a C(alpha)
        neighbors = [a.GetAtomicNum() for a in amide_n.GetNeighbors()]
        if 6 in neighbors:
            valid_connections +=1

    if valid_connections != 2:
        return False, f"Found {valid_connections} peptide bonds connected to amino acids, need exactly 2"

    # Check molecular features typical of peptides
    # Calculate the number of amide bonds
    num_amide_bonds = rdMolDescriptors.CalcNumAmideBonds(mol)
    if num_amide_bonds != 2:
        return False, f"Found {num_amide_bonds} amide bonds, expected 2 for a tripeptide"

    # Consider possible cyclization
    # If the molecule is cyclic, adjust the expected number of peptide bonds
    is_cyclic = mol.GetRingInfo().IsConnected()
    if is_cyclic and num_amide_bonds < 2:
        return False, "Molecule is cyclic but does not have enough peptide bonds"

    return True, "Molecule is a tripeptide with three amino acid residues connected by peptide bonds"