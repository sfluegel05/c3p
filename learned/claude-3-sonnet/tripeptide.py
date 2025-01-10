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
        return False, f"Found {len(peptide_matches)} peptide bonds, need exactly 2"
    
    # Look for terminal carboxyl group pattern (including salts/esters)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-,OX2]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Missing terminal carboxyl group or derivative"

    # Look for N-terminal pattern (includes both free NH2 and modified forms)
    nterm_pattern = Chem.MolFromSmarts("[$([NX3H2,NX3H3+]),$([NX3][CX3](=[OX1]))][CX4][CX3](=[OX1])[NX3]")
    if not mol.HasSubstructMatch(nterm_pattern):
        return False, "Missing N-terminal amino acid residue pattern"

    # Pattern for complete amino acid residue
    aa_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])")
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    
    # Should find at least 3 amino acid residues
    if len(aa_matches) < 3:
        return False, f"Found only {len(aa_matches)} amino acid residues"

    # Check for cyclic peptides (should be linear)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 2:  # Allow up to 2 rings for aromatic side chains
        return False, "Structure appears to be cyclic"

    # Count carbons, nitrogens and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 6:
        return False, f"Too few carbons ({c_count}) for tripeptide"
    if n_count < 3:
        return False, f"Too few nitrogens ({n_count}) for tripeptide"
    if o_count < 4:
        return False, f"Too few oxygens ({o_count}) for tripeptide"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:  # Minimum for smallest possible tripeptide
        return False, f"Molecular weight ({mol_wt:.1f}) too low for tripeptide"
    if mol_wt > 1000:  # Maximum for typical tripeptide
        return False, f"Molecular weight ({mol_wt:.1f}) too high for tripeptide"

    # Check for reasonable number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Too few rotatable bonds for tripeptide"
    
    # Additional check for backbone connectivity
    backbone_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=[OX1])[NX3][CX4][CX3](=[OX1])[NX3][CX4][CX3](=[OX1])")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Missing complete tripeptide backbone pattern"

    return True, "Contains characteristic tripeptide structure with 3 amino acid residues"