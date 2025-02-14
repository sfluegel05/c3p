"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to find a peptide bond
    peptide_bond_pattern = Chem.MolFromSmarts("[-N]-C(=O)-")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # SMARTS pattern to find amino acid residue
    # This pattern tries to capture amino acid residues
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]-[CX4;H2,H1,H0]-[CX3](=[OX1])")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    # If we find less than 10 amino acid residues, it cannot be a polypeptide
    if len(amino_acid_matches) < 10:
        return False, f"Molecule contains {len(amino_acid_matches)} amino acid residues, needs at least 10"
    
    # Check to see that the molecule contains at least 9 peptide bonds.
    # Although the pattern is not fool proof, the number of peptide bonds should be around 1 less than the number of amino acid residues.
    if len(peptide_bonds) < 9:
         return False, f"Molecule contains {len(peptide_bonds)} peptide bonds, which seems too low."


    return True, "Molecule contains at least 10 amino acid residues connected by peptide bonds"