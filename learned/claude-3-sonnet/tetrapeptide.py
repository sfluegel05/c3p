"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: tetrapeptide
A molecule containing exactly four amino acid residues connected by peptide linkages
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is a tetrapeptide and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate basic properties
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight > 1000:  # Most tetrapeptides are under 1000 Da
        return False, "Molecular weight too high for typical tetrapeptide"
    
    # Pattern for amino acid residues (including proline)
    aa_pattern = Chem.MolFromSmarts("[NX3,NX4;H0,H1,H2][CX4][CX3](=[OX1])[NX3,OX2H,OX1-,N]")
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    
    # Pattern for peptide bonds (including proline)
    peptide_pattern = Chem.MolFromSmarts("[$([NX3H1,NX4H2]),$([NX3](C)(C))][CX4][CX3](=[OX1])[NX3H1,NX4H2]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    # Pattern for N-terminus (including modified)
    n_term_pattern = Chem.MolFromSmarts("[$([NX3H2,NX4H3+]),$([NX3H1](C))][CX4][CX3](=O)")
    n_term = len(mol.GetSubstructMatches(n_term_pattern))
    
    # Pattern for C-terminus (including modified)
    c_term_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[$([OX2H]),$([OX1-]),$([NX3H2])]")
    c_term = len(mol.GetSubstructMatches(c_term_pattern))
    
    # Count amino acid residues
    aa_count = len(aa_matches)
    if aa_count != 4:
        return False, f"Found {aa_count} amino acid residues, need exactly 4"
    
    # Count peptide bonds
    peptide_count = len(peptide_matches)
    if peptide_count < 3:
        return False, f"Found {peptide_count} peptide bonds, need at least 3"
    
    # Check for cyclic structure
    ring_info = mol.GetRingInfo()
    ring_size = max([len(r) for r in ring_info.AtomRings()], default=0)
    
    # Exclude molecules with large rings (except for proline rings)
    if ring_size > 12:  # Proline ring is 5-membered
        return False, "Ring size too large for tetrapeptide"
        
    # Check for appropriate terminal groups (unless cyclic)
    if ring_size < 12:  # Linear peptide
        if n_term == 0 and c_term == 0:
            return False, "Missing both N-terminus and C-terminus"
    
    # Count backbone atoms to ensure proper chain length
    backbone_pattern = Chem.MolFromSmarts("[NX3,NX4][CX4][CX3]=O")
    backbone_atoms = len(mol.GetSubstructMatches(backbone_pattern))
    if backbone_atoms < 4:
        return False, "Insufficient backbone length for tetrapeptide"
        
    # Final classification
    if ring_size >= 12:
        return True, "Cyclic tetrapeptide with 4 amino acid residues connected by peptide bonds"
    else:
        return True, "Linear tetrapeptide with 4 amino acid residues connected by peptide bonds"