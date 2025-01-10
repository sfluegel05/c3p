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

    # Pattern for peptide bonds (including proline-type)
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])[NX3]")
    
    # More flexible pattern for amino acid residues
    aa_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])[*]")
    
    # Pattern for N-terminal amino group (including modifications)
    n_term_pattern = Chem.MolFromSmarts("[$([NX3H2]),$([NX3H](C=O)),$([NX3](C=O)C)][CX4][CX3]=O")
    
    # Pattern for C-terminal group (including modifications)
    c_term_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[$([OX2H]),$([OX2]C),$([NX3H2]),$([NX3H]C)]")
    
    # Get all matches
    peptide_bonds = mol.GetSubstructMatches(peptide_pattern)
    aa_residues = mol.GetSubstructMatches(aa_pattern)
    
    # Get ring info
    ring_info = mol.GetRingInfo()
    is_cyclic = ring_info.NumRings() > 0
    
    # Count unique amino acid residues (avoiding overlaps)
    unique_residues = set()
    for match in aa_residues:
        central_C = match[1]  # The alpha carbon
        if central_C not in unique_residues:
            unique_residues.add(central_C)
    
    aa_count = len(unique_residues)
    
    # Check sequence connectivity
    if not is_cyclic:
        # For linear peptides, check termini
        n_term_matches = mol.GetSubstructMatches(n_term_pattern)
        c_term_matches = mol.GetSubstructMatches(c_term_pattern)
        
        if not (n_term_matches and c_term_matches):
            return False, "Missing required N- or C-terminus"
        
        # Check for continuous peptide backbone
        backbone = set()
        for bond in peptide_bonds:
            backbone.update(bond)
        
        if len(backbone) < 12:  # Minimum atoms for tetrapeptide backbone
            return False, "Insufficient peptide backbone connectivity"
    else:
        # For cyclic peptides, check ring size and connectivity
        rings = ring_info.AtomRings()
        max_ring_size = max(len(ring) for ring in rings)
        
        if max_ring_size < 12 or max_ring_size > 16:
            return False, "Ring size not consistent with tetrapeptide"
    
    # Count peptide bonds (avoiding overlaps)
    unique_peptide_bonds = set()
    for bond in peptide_bonds:
        if bond[0] not in unique_peptide_bonds and bond[-1] not in unique_peptide_bonds:
            unique_peptide_bonds.update(bond)
    
    peptide_bond_count = len(peptide_bonds)
    
    # Classification logic
    if aa_count != 4:
        return False, f"Found {aa_count} amino acid residues, need exactly 4"
    
    if peptide_bond_count < 3:
        return False, f"Found only {peptide_bond_count} peptide bonds, need at least 3"
    
    if peptide_bond_count > 4:
        return False, f"Found {peptide_bond_count} peptide bonds, too many for tetrapeptide"
    
    # Success cases
    if is_cyclic:
        return True, "Cyclic tetrapeptide with 4 amino acid residues"
    else:
        return True, "Linear tetrapeptide with 4 amino acid residues"