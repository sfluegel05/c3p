"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: alpha-amino acid
An amino acid in which the amino group is located on the carbon atom at the position alpha to the carboxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens for better pattern matching
    mol = Chem.AddHs(mol)
    
    # Look for carboxylic acid group pattern
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Find all carboxyl groups
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # Look for carbon with primary amine (-NH2) pattern
    # [CX4] carbon with 4 bonds
    # [NX3H2] nitrogen with 3 bonds and 2 hydrogens (primary amine)
    alpha_amino_pattern = Chem.MolFromSmarts("[NX3H2][CX4]C(=O)[OH]")
    
    if not mol.HasSubstructMatch(alpha_amino_pattern):
        return False, "No alpha amino acid pattern found"
    
    # Additional check for substituted amines (still alpha amino acids)
    alpha_amino_pattern2 = Chem.MolFromSmarts("[NX3H1;!$(NC=O)][CX4]C(=O)[OH]")
    alpha_amino_pattern3 = Chem.MolFromSmarts("[NX3H0;!$(NC=O)][CX4]C(=O)[OH]")
    
    # Get all matches
    matches1 = mol.GetSubstructMatches(alpha_amino_pattern)
    matches2 = mol.GetSubstructMatches(alpha_amino_pattern2)
    matches3 = mol.GetSubstructMatches(alpha_amino_pattern3)
    
    total_matches = len(matches1) + len(matches2) + len(matches3)
    
    if total_matches == 0:
        return False, "No alpha carbon with amino group found"
    
    # Check if the amino group is not part of an amide
    amide_pattern = Chem.MolFromSmarts("[NX3H2]C(=O)")
    peptide_pattern = Chem.MolFromSmarts("[NX3H1]C(=O)")
    
    # If all nitrogens are part of amides/peptides, it's not a free alpha amino acid
    if mol.HasSubstructMatch(amide_pattern) or mol.HasSubstructMatch(peptide_pattern):
        # Need to verify that there's at least one non-amide amino group
        non_amide_amino = False
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7:  # Nitrogen
                if atom.GetDegree() <= 3 and not any(n.GetAtomicNum() == 6 and any(b.GetBondType() == Chem.BondType.DOUBLE for b in n.GetBonds()) for n in atom.GetNeighbors()):
                    non_amide_amino = True
                    break
        if not non_amide_amino:
            return False, "All amino groups are part of amides/peptides"
    
    # Remove explicit hydrogens for final checks
    mol = Chem.RemoveHs(mol)
    
    # Count basic elements
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 2:
        return False, "Too few carbons for an alpha amino acid"
    if n_count < 1:
        return False, "No nitrogen found"
    if o_count < 2:
        return False, "Too few oxygens for carboxylic acid group"
        
    return True, "Contains alpha carbon with amino group and carboxylic acid"