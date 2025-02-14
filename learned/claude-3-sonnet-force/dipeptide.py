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
    
    # Look for peptide bond pattern (-CO-NH-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)NCC(=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 1:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, expected exactly 1"
    
    # Look for two amino acid residues (atoms with both -NH2 and -COOH groups)
    residue_pattern = Chem.MolFromSmarts("NC(C(=O)O)C")
    residue_matches = mol.GetSubstructMatches(residue_pattern)
    if len(residue_matches) != 2:
        return False, f"Found {len(residue_matches)} amino acid residues, expected exactly 2"
    
    # Ensure the peptide bond connects the two residues
    peptide_bond_atoms = set(atom.GetIdx() for atom in mol.GetAtoms()[peptide_bond_matches[0]])
    residue_atoms = set(sum([list(mol.GetAtoms()[res_match]) for res_match in residue_matches], []))
    if not peptide_bond_atoms.issubset(residue_atoms):
        return False, "Peptide bond does not connect the two amino acid residues"
    
    return True, "Contains two amino acid residues connected by a peptide bond"