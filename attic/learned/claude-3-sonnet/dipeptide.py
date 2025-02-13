"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide contains two amino acid residues connected by a peptide linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible peptide bond pattern
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    if len(peptide_matches) == 0:
        return False, "No peptide bond found"
        
    if len(peptide_matches) > 3:
        return False, "Too many peptide bonds for dipeptide"

    # Check for cyclic dipeptides (diketopiperazines)
    diketopiperazine = Chem.MolFromSmarts("O=C1NC(=O)CN1")
    if mol.HasSubstructMatch(diketopiperazine):
        return True, "Cyclic dipeptide (diketopiperazine) structure found"
    
    # More flexible patterns for terminal groups
    terminal_amine_patterns = [
        "[NX3H2,NX4H3+][CX4]",  # Primary amine
        "[NX3H1,NX4H2+][CX4]",  # Secondary amine
        "[NX3H0,NX4H1+][CX4]"   # Tertiary amine
    ]
    
    terminal_carboxyl_patterns = [
        "[CX3](=[OX1])[OX2H1,OX1-]",  # Free acid
        "[CX3](=[OX1])[OX2]",         # Ester
        "[CX3](=[OX1])[NX3]"          # Amide
    ]
    
    # Check for any terminal group patterns
    has_terminal_amine = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) 
                           for pattern in terminal_amine_patterns)
    has_terminal_carboxyl = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) 
                              for pattern in terminal_carboxyl_patterns)

    # Alpha carbon pattern (more flexible)
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])")
    alpha_carbons = mol.GetSubstructMatches(alpha_carbon_pattern)
    
    # Basic requirements for linear dipeptide
    if len(peptide_matches) >= 1 and len(peptide_matches) <= 2:
        if len(alpha_carbons) >= 1:
            # Count C, N, O atoms
            c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
            o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            
            if c_count < 3:
                return False, "Too few carbons for dipeptide"
            if n_count < 1:
                return False, "Too few nitrogens for dipeptide"
            if o_count < 2:
                return False, "Too few oxygens for dipeptide"
            
            # More flexible molecular weight range
            mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
            if mol_wt < 80:  # Relaxed lower bound
                return False, "Molecular weight too low for dipeptide"
                
            # If we have a peptide bond and reasonable atom counts, likely a dipeptide
            if has_terminal_amine or has_terminal_carboxyl:
                return True, "Contains peptide bond with appropriate terminal groups"
            else:
                # Still might be a modified dipeptide
                return True, "Contains peptide bond with modified terminal groups"
            
    return False, "Does not match dipeptide structure requirements"