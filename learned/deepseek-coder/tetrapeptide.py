"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is a molecule that contains four amino-acid residues connected by peptide linkages.

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

    # More inclusive amino acid residue pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H0,H1,H2][CX4H](*)[CX3](=[OX1])[OX2H0,H1]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # A tetrapeptide should have exactly 4 amino acid residues
    if len(amino_acid_matches) != 4:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, need exactly 4"

    # More specific peptide bond pattern that connects amino acid residues
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0,H1][CX4H](*)[CX3](=[OX1])")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # A tetrapeptide should have exactly 3 peptide bonds
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3"

    # Check molecular weight - tetrapeptides typically have a molecular weight between 400-800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for tetrapeptides"

    # Check for linearity of the peptide chain
    # Count the number of rotatable bonds - should be between 3-10 for a tetrapeptide
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3 or n_rotatable > 10:
        return False, f"Number of rotatable bonds ({n_rotatable}) is outside the typical range for tetrapeptides"

    # Additional check for terminal groups
    has_n_terminal = mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3H2][CX4H](*)[CX3](=[OX1])"))
    has_c_terminal = mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3H0,H1][CX4H](*)[CX3](=[OX1])[OX2H1]"))
    if not (has_n_terminal and has_c_terminal):
        return False, "Missing proper N-terminal or C-terminal groups"

    return True, "Contains four amino acid residues connected by exactly three peptide bonds with proper terminal groups"