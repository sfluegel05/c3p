"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: oligopeptide
A peptide containing a relatively small number of amino acids.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for peptide bonds (-C(=O)-NH-)
    peptide_pattern = Chem.MolFromSmarts("[NX3H1][CX3](=[OX1])[CX4]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    if len(peptide_matches) == 0:
        return False, "No peptide bonds found"
        
    # Count number of peptide bonds to estimate number of amino acids
    num_peptide_bonds = len(peptide_matches)
    
    # Look for characteristic amino acid patterns
    amino_pattern = Chem.MolFromSmarts("[NX3H2,NX3H1][CX4H1][CX3](=[OX1])")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    if len(amino_matches) == 0:
        return False, "No amino acid residues found"
        
    # Count atoms to estimate size
    num_atoms = mol.GetNumAtoms()
    if num_atoms > 300:  # Rough upper limit for oligopeptides
        return False, "Molecule too large for oligopeptide"
        
    # Count nitrogens and oxygens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 2 or o_count < 2:
        return False, "Insufficient N/O atoms for peptide structure"
    
    # Estimate number of amino acids (rough approximation)
    estimated_aas = max(num_peptide_bonds + 1, len(amino_matches))
    
    if estimated_aas > 50:
        return False, f"Too many amino acids ({estimated_aas}) for oligopeptide"
    elif estimated_aas < 2:
        return False, "Too few amino acids for oligopeptide"
        
    # Look for carboxyl group (may be modified in some cases)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1,OX1-,OX2]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        # Check if it's a cyclic peptide
        ring_info = mol.GetRingInfo()
        if not ring_info.NumRings():
            return False, "No carboxyl group or cyclic structure found"
    
    # Additional check for characteristic peptide molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:  # Smaller than dipeptide
        return False, "Molecular weight too low for oligopeptide"
        
    # If we've made it here, it's likely an oligopeptide
    return True, f"Contains {estimated_aas} amino acids connected by peptide bonds"