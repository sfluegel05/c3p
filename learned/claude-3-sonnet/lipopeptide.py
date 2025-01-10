"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bonds with more inclusive patterns
    peptide_patterns = [
        "[NX3][CX3](=[OX1])[CX4]",  # Standard peptide bond
        "[NX3;H0,H1][CX3](=[OX1])[CX4]",  # Include N-methylated
        "[NX3][CX3](=[OX1])[CX4,CX3]"  # More flexible C-terminal
    ]
    
    peptide_count = 0
    for pattern in peptide_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            peptide_count += len(mol.GetSubstructMatches(pat))
    
    if peptide_count < 2:
        return False, "Insufficient peptide bonds found"

    # Look for amino acid residues with more patterns
    aa_patterns = [
        "[NX3H1,NX3H2][CX4H]([CX4,CX3])[CX3](=[OX1])[OX2,NX3]",  # Standard
        "[NX3][CX4H]([CX4,CX3])[CX3](=[OX1])[OX2,NX3]",  # Modified
        "[NX3][CX4H]([*])[CX3](=[OX1])[OX2,NX3]",  # General
        "[NX3]([CH3])[CX4H]([*])[CX3](=[OX1])[OX2,NX3]"  # N-methylated
    ]
    
    aa_count = 0
    for pattern in aa_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            aa_count += len(mol.GetSubstructMatches(pat))
    
    if aa_count < 2:
        return False, "Insufficient amino acid residues found"

    # Enhanced lipid chain detection
    lipid_patterns = [
        "[CX4H2][CX4H2][CX4H2][CX4H2][CX4H2]",  # Shorter chain (C5)
        "[CX4H2,CX4H][CX4H2,CX4H][CX4H2,CX4H][CX4H2,CX4H]",  # Branched
        "[CH3][CX4H2][CX4H2][CX4H2][CX4H2]",  # Terminal methyl
        "[$([CX4H2][CX4H2][CX4H2][CX4H2][F,Cl,Br,I])]",  # Halogenated
        "[$([CX4H2][CX4H2][CX4H2][CX4H2][CF,CCl,CBr,CI])]"  # Halogenated branch
    ]
    
    has_lipid = False
    lipid_size = 0
    for pattern in lipid_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            matches = mol.GetSubstructMatches(pat)
            if matches:
                has_lipid = True
                lipid_size = max(lipid_size, len(matches[0]) if matches else 0)
    
    if not has_lipid or lipid_size < 5:  # Reduced minimum size
        return False, "No suitable lipid chain found"

    # Count atoms and check ratios
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 2 or o_count < 2:
        return False, "Insufficient N/O atoms for peptide structure"
    
    if c_count < 10:  # Reduced minimum
        return False, "Carbon count too low for lipopeptide"

    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if mol_wt < 350:  # Reduced threshold
        return False, "Molecular weight too low for lipopeptide"
    
    if rotatable_bonds < 8:  # Reduced threshold
        return False, "Insufficient flexibility for lipopeptide"

    # Improved cyclic peptide detection
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Look for peptide cycles with stricter requirements
    has_peptide_ring = False
    if any(size >= 6 for size in ring_sizes):
        cycle_patterns = [
            "[NX3][CX3](=[OX1])[CX4]~[NX3][CX3](=[OX1])[CX4]~[NX3]",  # 3+ peptide bonds
            "[NX3][CX3](=[OX1])[CX4]~[NX3][CX3](=[OX1])[CX4]~[NX3][CX3](=[OX1])"  # Alternative
        ]
        for pattern in cycle_patterns:
            pat = Chem.MolFromSmarts(pattern)
            if pat is not None and mol.HasSubstructMatch(pat):
                has_peptide_ring = True
                break
    
    # Check peptide bond to ring size ratio for cyclic structures
    if has_peptide_ring:
        max_ring_size = max(ring_sizes)
        if peptide_count >= max_ring_size/3:  # At least 1/3 of ring bonds should be peptide bonds
            return True, "Cyclic lipopeptide structure identified"
    
    # Check for branched peptide structure
    branch_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]([NX3])[CX3](=[OX1])")
    if branch_pattern is not None and mol.HasSubstructMatch(branch_pattern):
        return True, "Branched lipopeptide structure identified"

    # Linear peptide validation
    if peptide_count >= 3 and aa_count >= 3:  # Require more evidence for linear peptides
        return True, "Linear lipopeptide with sufficient peptide bonds and lipid chain"
        
    return False, "Does not meet lipopeptide structural requirements"