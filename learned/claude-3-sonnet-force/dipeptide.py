"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: CHEBI:36556 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule containing two amino-acid residues connected by peptide linkages.

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
    
    # Look for amino acid residues (atoms with both -NH2 and -COOH groups)
    amino_acid_pattern = Chem.MolFromSmarts("[NH2,NH3][CX4H][CX3](=[OX1])[OX2H,OX1-]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) != 2:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, expected exactly 2"
    
    # Look for peptide bond patterns (-CO-NH- or -CO-N[C@@H]-)
    peptide_bond_pattern_1 = Chem.MolFromSmarts("C(=O)NC")
    peptide_bond_pattern_2 = Chem.MolFromSmarts("C(=O)N[C@@H]")
    peptide_bond_matches_1 = mol.GetSubstructMatches(peptide_bond_pattern_1)
    peptide_bond_matches_2 = mol.GetSubstructMatches(peptide_bond_pattern_2)
    if len(peptide_bond_matches_1) + len(peptide_bond_matches_2) != 1:
        return False, f"Found {len(peptide_bond_matches_1) + len(peptide_bond_matches_2)} peptide bonds, expected exactly 1"
    
    # Ensure the peptide bond connects the two amino acid residues
    peptide_bond_atoms = set()
    if peptide_bond_matches_1:
        peptide_bond_atoms = set(atom.GetIdx() for atom in mol.GetAtoms()[peptide_bond_matches_1[0]])
    else:
        peptide_bond_atoms = set(atom.GetIdx() for atom in mol.GetAtoms()[peptide_bond_matches_2[0]])
    
    residue_atoms = set(sum([list(mol.GetAtoms()[match]) for match in amino_acid_matches], []))
    if not peptide_bond_atoms.issubset(residue_atoms):
        return False, "Peptide bond does not connect the two amino acid residues"
    
    return True, "Contains two amino acid residues connected by a peptide bond"