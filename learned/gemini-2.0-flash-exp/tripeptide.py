"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is a peptide consisting of three amino acid residues connected by two peptide bonds.

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
    
    # Check for two peptide bonds - modified to look for non-terminal peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX2][CX4]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 2:
         return False, f"Found {len(peptide_bond_matches)} peptide bonds, expected 2"

    # Amino acid residue pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX2][CX4]([H])C(=O)[CX3]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) != 3:
         return False, f"Found {len(amino_acid_matches)} amino acid residues, expected 3."
    
    # check for free carboxylic acid
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) < 1:
         return False, f"Found {len(carboxylic_acid_matches)} free carboxyl groups, expected at least 1"

    # Count number of carbon atoms: must be in the 10-25 range for the size of most tripeptides.
    num_carbon_atoms = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if not 10 <= num_carbon_atoms <= 30:
        return False, f"Number of carbon atoms {num_carbon_atoms} not in expected range (10-30)"
    
    # Count number of oxygen atoms: must be in the range of 4-8 for most tripeptides.
    num_oxygen_atoms = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if not 4 <= num_oxygen_atoms <= 9:
        return False, f"Number of oxygen atoms {num_oxygen_atoms} not in expected range (4-9)"
    
    # check for molecular weight (should be between 200 and 600)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if not 200 <= mol_weight <= 700:
         return False, f"Molecular weight {mol_weight} not in expected range (200-700)"

    return True, "Contains three amino acid residues connected by two peptide bonds."