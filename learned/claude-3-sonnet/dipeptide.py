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

    # Look for peptide bonds (-C(=O)-NH-)
    peptide_pattern = Chem.MolFromSmarts("[NX3H1][CX3](=[OX1])[CX4]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    if len(peptide_matches) < 1:
        return False, "No peptide bond found"
    
    # For cyclic dipeptides, look for diketopiperazine pattern
    diketopiperazine = Chem.MolFromSmarts("O=C1NC(=O)CN1")
    if mol.HasSubstructMatch(diketopiperazine):
        return True, "Cyclic dipeptide (diketopiperazine) structure found"
    
    # For linear dipeptides:
    # Look for terminal amine (-NH2) and carboxyl (-COOH) groups
    terminal_amine = Chem.MolFromSmarts("[NX3H2,NX4H3+][CX4]")
    terminal_carboxyl = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1,OX1-]")
    
    amine_matches = mol.GetSubstructMatches(terminal_amine)
    carboxyl_matches = mol.GetSubstructMatches(terminal_carboxyl)
    
    # Count alpha carbons (carbons next to peptide bond nitrogen)
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3H1][CX4][CX3](=[OX1])")
    alpha_carbons = mol.GetSubstructMatches(alpha_carbon_pattern)
    
    # Basic requirements for linear dipeptide
    if len(peptide_matches) == 1 and len(alpha_carbons) >= 1:
        if len(amine_matches) >= 1 and len(carboxyl_matches) >= 1:
            # Check molecular weight - most dipeptides are between 150-400 Da
            mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
            if mol_wt < 100:
                return False, "Molecular weight too low for dipeptide"
            
            # Count C, N, O atoms
            c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
            o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            
            if c_count < 4:
                return False, "Too few carbons for dipeptide"
            if n_count < 2:
                return False, "Too few nitrogens for dipeptide"
            if o_count < 3:
                return False, "Too few oxygens for dipeptide"
                
            return True, "Contains peptide bond with terminal amine and carboxyl groups"
            
    # If more than one peptide bond, likely a larger peptide
    if len(peptide_matches) > 2:
        return False, "Too many peptide bonds for dipeptide"
        
    return False, "Does not match dipeptide structure requirements"