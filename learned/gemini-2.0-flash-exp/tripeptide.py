"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is a peptide consisting of three amino acid residues connected by peptide bonds.

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
    
    # Check for two peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 2:
         return False, f"Found {len(peptide_bond_matches)} peptide bonds, expected 2"

    # Amino acid pattern (simple, requires an alpha-C with an amine and a carbonyl)
    amino_acid_pattern = Chem.MolFromSmarts("[CX4]([NX2])C(=O)[OX1]") # Modified to reflect -C(=O)O
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    # Check for three amino acid residues
    if len(amino_acid_matches) < 3:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, expected 3."

    # Count number of carbon atoms: must be in the 10-25 range for the size of most tripeptides.
    num_carbon_atoms = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if not 10 <= num_carbon_atoms <= 25:
        return False, f"Number of carbon atoms {num_carbon_atoms} not in expected range (10-25)"
    
    # Count number of oxygen atoms: must be in the range of 4-7 for most tripeptides.
    num_oxygen_atoms = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if not 4 <= num_oxygen_atoms <= 7:
        return False, f"Number of oxygen atoms {num_oxygen_atoms} not in expected range (4-7)"
    
    # check for molecular weight (should be between 200 and 600)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if not 200 <= mol_weight <= 600:
         return False, f"Molecular weight {mol_weight} not in expected range (200-600)"

    return True, "Contains three amino acid residues connected by two peptide bonds."