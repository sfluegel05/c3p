"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: tetrapeptide
A molecule containing exactly four amino acid residues connected by peptide linkages
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Look for peptide bonds (-C(=O)-NH-)
    peptide_pattern = Chem.MolFromSmarts("[NX3,NX4;H1,H2][CX4][CX3](=[OX1])[NX3,NX4;H1,H2]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    # For tetrapeptide, we expect 3 peptide bonds (connecting 4 amino acids)
    if len(peptide_matches) < 3:
        return False, f"Found only {len(peptide_matches)} peptide bonds, need 3 for tetrapeptide"
        
    # Look for alpha carbons with attached NH group
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3,NX4;H1,H2][CX4][CX3](=O)")
    alpha_carbons = mol.GetSubstructMatches(alpha_carbon_pattern)
    
    # We expect 4 alpha carbons for a tetrapeptide
    if len(alpha_carbons) < 4:
        return False, f"Found only {len(alpha_carbons)} alpha carbons, need 4 for tetrapeptide"
    
    # Look for N-terminus (primary or secondary amine)
    n_terminus_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3+][CX4][CX3](=O)")
    n_terminus = mol.GetSubstructMatches(n_terminus_pattern)
    
    # Look for C-terminus (carboxyl or similar group)
    c_terminus_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-,N]")
    c_terminus = mol.GetSubstructMatches(c_terminus_pattern)
    
    # Check for minimum required terminal groups
    if not (n_terminus or c_terminus):
        return False, "Missing both N-terminus and C-terminus"
        
    # Count nitrogens that are part of peptide bonds or terminal groups
    peptide_n_pattern = Chem.MolFromSmarts("[NX3,NX4;H1,H2][CX4][CX3]=O")
    peptide_nitrogens = len(mol.GetSubstructMatches(peptide_n_pattern))
    
    # For a tetrapeptide, we expect 4 nitrogens involved in peptide bonds/terminals
    if peptide_nitrogens < 4:
        return False, f"Found only {peptide_nitrogens} peptide nitrogens, need 4 for tetrapeptide"
    
    # Count carbonyls that are part of peptide bonds or terminal groups
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3,OX2H,OX1-]")
    carbonyls = len(mol.GetSubstructMatches(carbonyl_pattern))
    
    # For a tetrapeptide, we expect 4 carbonyls (3 peptide bonds + C-terminus)
    if carbonyls < 4:
        return False, f"Found only {carbonyls} peptide carbonyls, need 4 for tetrapeptide"
    
    # Additional check for cyclic peptides
    cycle_pattern = Chem.MolFromSmarts("[NX3,NX4;H1,H2][CX4][CX3](=[OX1])[NX3,NX4;H1,H2]")
    cycle_matches = mol.GetSubstructMatches(cycle_pattern)
    is_cyclic = len(cycle_matches) >= 4  # If we find 4 or more peptide bonds, it might be cyclic
    
    # Final classification
    if is_cyclic:
        return True, "Cyclic tetrapeptide with 4 amino acid residues connected by peptide bonds"
    else:
        return True, "Linear tetrapeptide with 4 amino acid residues connected by peptide bonds"