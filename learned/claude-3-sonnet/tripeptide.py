"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:47926 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of exactly 3 amino acid residues connected by peptide bonds.

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

    # Basic peptide bond pattern (amide bond)
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    if peptide_pattern is None:
        return False, "Invalid peptide SMARTS pattern"
    
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    # For a tripeptide, we need exactly 2 peptide bonds
    if len(peptide_matches) < 2:
        return False, f"Found {len(peptide_matches)} peptide bonds, need at least 2"

    # Look for amino group pattern
    amino_pattern = Chem.MolFromSmarts("[NX3H2,NX3H3+][CX4]")
    if amino_pattern is None:
        return False, "Invalid amino SMARTS pattern"
        
    # Look for carboxyl group pattern
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-]")
    if carboxyl_pattern is None:
        return False, "Invalid carboxyl SMARTS pattern"
    
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "Missing terminal amino group"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Missing terminal carboxyl group"

    # Count nitrogens and oxygens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Minimum requirements for tripeptide:
    # - At least 3 nitrogens (1 terminal + 2 peptide bonds)
    # - At least 4 oxygens (2 peptide C=O + 1 terminal COOH)
    if n_count < 3:
        return False, f"Too few nitrogens ({n_count}) for a tripeptide"
    if o_count < 4:
        return False, f"Too few oxygens ({o_count}) for a tripeptide"

    # Check molecular weight - more permissive range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:  # More permissive minimum weight
        return False, f"Molecular weight ({mol_wt:.1f}) too low for tripeptide"
    if mol_wt > 2000:  # More permissive maximum weight
        return False, f"Molecular weight ({mol_wt:.1f}) too high for tripeptide"

    # Look for alpha-amino acid pattern (including proline-like residues)
    aa_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=[OX1])")
    if aa_pattern is None:
        return False, "Invalid amino acid SMARTS pattern"
        
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    
    # Should find at least 2 internal amino acid residues
    # (the terminal one might not match this pattern exactly)
    if len(aa_matches) < 2:
        return False, f"Found only {len(aa_matches)} amino acid residues"

    # Check for reasonable number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:  # Minimum expected for tripeptide
        return False, "Too few rotatable bonds for tripeptide"

    # Count carbon atoms (should have at least 6 for smallest possible tripeptide)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, f"Too few carbons ({c_count}) for tripeptide"

    return True, "Contains characteristic tripeptide structure with 3 amino acid residues"