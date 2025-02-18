"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:36357 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import Bond

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide consisting of three amino acid residues
    connected by peptide bonds.

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

    # Count peptide bonds (-C(=O)-N-)
    peptide_bonds = sum(
        1 for bond in mol.GetBonds() if bond.GetBondType() == Bond.AMIDE
    )
    if peptide_bonds != 2:
        return False, f"Found {peptide_bonds} peptide bonds, need exactly 2"

    # Count amino acid residues (-NH-C(=O)-C-)
    amino_acids = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "N" and atom.GetHybridization() == Chem.HybridizationType.SP3:
            if atom.GetBondRingSize(Bond.AMIDE) == 0:  # Ignore ring systems
                amino_acids += 1
    if amino_acids != 3:
        return False, f"Found {amino_acids} amino acid residues, need exactly 3"

    # Check for atoms other than C, N, O, H, S
    allowed_atoms = {"C", "N", "O", "H", "S"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowed_atoms:
            return False, "Contains atoms other than C, N, O, H, S"

    return True, "Contains three amino acid residues connected by peptide bonds"