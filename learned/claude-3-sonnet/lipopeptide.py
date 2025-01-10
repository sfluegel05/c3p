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

    # Look for peptide bonds (-C(=O)-N-) with connected carbons
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    if peptide_pattern is None:
        return None, "Invalid peptide SMARTS pattern"
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    if len(peptide_matches) < 2:  # Need at least 2 peptide bonds
        return False, "Insufficient peptide bonds found"

    # Look for amino acid residues
    aa_patterns = [
        "[NX3H1,NX3H2][CX4H]([CX4,CX3])[CX3](=[OX1])[OX2,NX3]",  # Standard amino acid
        "[NX3][CX4H]([CX4,CX3])[CX3](=[OX1])[OX2,NX3]",          # Modified amino acid
        "[NX3][CX4H]([*])[CX3](=[OX1])[OX2,NX3]"                 # General amino acid
    ]
    
    aa_count = 0
    for pattern in aa_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            aa_count += len(mol.GetSubstructMatches(pat))
    
    if aa_count < 2:
        return False, "Insufficient amino acid residues found"

    # Look for lipid chains - using multiple patterns
    lipid_patterns = [
        "[CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2]",  # Straight chain
        "[CX4H2,CX4H][CX4H2,CX4H][CX4H2,CX4H][CX4H2,CX4H][CX4H2,CX4H]",  # Branched
        "[CH3][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2]"  # Terminal methyl
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
    
    if not has_lipid or lipid_size < 6:
        return False, "No suitable lipid chain found"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 2:
        return False, "Insufficient nitrogen atoms for peptide structure"
    
    if c_count < 12:
        return False, "Carbon count too low for lipopeptide"

    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if mol_wt < 400:
        return False, "Molecular weight too low for lipopeptide"
    
    if rotatable_bonds < 10:
        return False, "Insufficient flexibility for lipopeptide"

    # Check for cyclic peptide structure
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Look for peptide cycles (typically 6+ atoms)
    has_large_ring = any(size >= 6 for size in ring_sizes)
    
    if has_large_ring:
        # Verify ring contains peptide bonds
        cycle_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]~[NX3][CX3](=[OX1])[CX4]")
        if cycle_pattern is not None and mol.HasSubstructMatch(cycle_pattern):
            return True, "Cyclic lipopeptide structure identified"
    
    # Check for branched peptide structure
    branch_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]([NX3])[CX3](=[OX1])")
    if branch_pattern is not None and mol.HasSubstructMatch(branch_pattern):
        return True, "Branched lipopeptide structure identified"

    # If we've made it here, we have a linear peptide with lipid
    return True, "Linear lipopeptide with sufficient peptide bonds and lipid chain"