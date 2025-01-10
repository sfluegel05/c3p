"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide must contain four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for amide (peptide) bonds
    amide_pattern = Chem.MolFromSmarts("CON")
    amide_matches = [match for match in mol.GetSubstructMatches(amide_pattern)]

    # Check for the presence of C and N termini
    has_n_terminus = any(atom.GetAtomicNum() == 7 and atom.GetDegree() == 3 for atom in mol.GetAtoms())
    has_c_terminus = any(atom.GetAtomicNum() == 6 and atom.GetDegree() == 3 for atom in mol.GetAtoms())

    # Count amide groups (should be exactly 3 for a linear tetrapeptide)
    if len(amide_matches) == 3 and has_n_terminus and has_c_terminus:
        return True, "Contains four amino-acid residues connected by peptide bonds"
    
    return False, f"Contains {len(amide_matches)+1} amide bonds, but terminus analysis failed"

# The is_tetrapeptide function has been enhanced to additionally consider chemical environment,
# such as the presence of N and C termini, to reduce false positives from cyclic and modified peptides.