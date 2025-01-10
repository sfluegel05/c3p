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

    # More inclusive pattern for peptide bonds
    peptide_pattern = Chem.MolFromSmarts("[NX3;H1,H2][CX4][CX3](=[OX1])[NX3;H1,H2]")
    
    # Pattern for proline-type peptide bonds
    proline_pattern = Chem.MolFromSmarts("[NX3;R][CX4][CX3](=[OX1])[NX3;H1,H2]")
    
    # Pattern for amino acid residues (more specific)
    aa_pattern = Chem.MolFromSmarts("[NX3;H1,H2,H3][CX4][CX3](=[OX1])[OX2,NX3]")
    
    # Get matches
    peptide_bonds = mol.GetSubstructMatches(peptide_pattern)
    proline_bonds = mol.GetSubstructMatches(proline_pattern)
    aa_residues = mol.GetSubstructMatches(aa_pattern)
    
    # Total peptide bonds (including proline)
    total_peptide_bonds = len(peptide_bonds) + len(proline_bonds)
    
    # Check ring info
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Identify peptide rings (containing peptide bonds)
    peptide_rings = []
    for ring in rings:
        ring_atoms = set(ring)
        for bond in peptide_bonds:
            if set(bond).issubset(ring_atoms):
                peptide_rings.append(ring)
                break
    
    is_cyclic = len(peptide_rings) > 0
    
    # Count amino acid residues
    aa_count = len(aa_residues)
    
    # Additional check for N and C termini in linear peptides
    n_term_pattern = Chem.MolFromSmarts("[NX3;H1,H2][CX4][CX3]=O")
    c_term_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2,NX3]")
    
    has_n_term = len(mol.GetSubstructMatches(n_term_pattern)) > 0
    has_c_term = len(mol.GetSubstructMatches(c_term_pattern)) > 0
    
    # Classification logic
    if aa_count < 4:
        return False, f"Found only {aa_count} amino acid residues, need 4"
    
    if total_peptide_bonds < 3:
        return False, f"Found only {total_peptide_bonds} peptide bonds, need at least 3"
    
    if is_cyclic:
        # Check if the peptide ring has appropriate size for tetrapeptide
        peptide_ring_size = max(len(ring) for ring in peptide_rings)
        if peptide_ring_size < 12 or peptide_ring_size > 16:
            return False, "Peptide ring size not consistent with tetrapeptide"
            
        # Verify the ring contains peptide bonds
        peptide_bonds_in_ring = sum(1 for ring in peptide_rings 
                                  for bond in peptide_bonds 
                                  if set(bond).issubset(set(ring)))
        if peptide_bonds_in_ring < 3:
            return False, "Insufficient peptide bonds in ring"
            
        return True, "Cyclic tetrapeptide with 4 amino acid residues"
    else:
        # Linear peptide checks
        if not (has_n_term and has_c_term):
            return False, "Missing required N- or C-terminus"
            
        # Check for sequential connectivity
        backbone_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)[NX3]")
        backbone_count = len(mol.GetSubstructMatches(backbone_pattern))
        if backbone_count < 3:
            return False, "Insufficient peptide backbone connectivity"
        
        return True, "Linear tetrapeptide with 4 amino acid residues"