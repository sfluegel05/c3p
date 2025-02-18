"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with an attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify peptide bonds (amide linkages between amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")  # Pattern for an amide bond
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) < 2:
        return False, "Insufficient peptide bonds found"

    # Identify amino acid residues (alpha-amino acids)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4H]([#6])[CX3](=O)")  # N-C-C(=O) pattern
    amino_acids = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acids) < 2:
        return False, "Insufficient amino acid residues found"

    # Identify lipid chains (long aliphatic chains)
    lipid_chain_pattern = Chem.MolFromSmarts("C[C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]")  # Chain of at least 7 carbons
    lipid_chains = mol.GetSubstructMatches(lipid_chain_pattern)
    if len(lipid_chains) == 0:
        return False, "No lipid chains found"

    # Check for sufficient length of lipid chain
    aliphatic_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 2 and not atom.IsInRing())
    if aliphatic_carbons < 7:
        return False, "Lipid chain is too short"

    # Check if both peptide and lipid parts are connected
    if not Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False):
        return False, "Molecule is disconnected"

    return True, "Contains both peptide bonds and lipid chains indicative of a lipopeptide"