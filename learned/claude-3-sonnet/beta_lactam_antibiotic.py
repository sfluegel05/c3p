"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: beta-lactam antibiotics
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    Beta-lactam antibiotics contain a beta-lactam ring and are used as antibiotics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for beta-lactam ring (4-membered ring with N and C=O)
    beta_lactam_pattern = Chem.MolFromSmarts("[N]1[C]([C][C]1)=O")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Count beta-lactam rings
    beta_lactam_matches = len(mol.GetSubstructMatches(beta_lactam_pattern))
    if beta_lactam_matches > 1:
        return False, "Multiple beta-lactam rings found - unusual for antibiotics"

    # Look for characteristic groups often found in beta-lactam antibiotics
    
    # Carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for common ring systems in beta-lactams
    
    # Penicillin-like thiazolidine ring fused to beta-lactam
    penam_pattern = Chem.MolFromSmarts("[N]1[C]2[C](S[C]([C]2)([C])[C])[C]1=O")
    
    # Cephalosporin-like 6-membered ring fused to beta-lactam
    cephem_pattern = Chem.MolFromSmarts("[N]1[C]2[C](S[C]=[C]2)[C]1=O")
    
    # Carbapenem-like 5-membered ring fused to beta-lactam
    carbapenem_pattern = Chem.MolFromSmarts("[N]1[C]2[C](=[C][C]2)[C]1=O")
    
    # Check for at least one of the common ring systems
    has_common_ring = any(
        mol.HasSubstructMatch(pattern) 
        for pattern in [penam_pattern, cephem_pattern, carbapenem_pattern]
    )
    
    if not has_common_ring:
        return False, "Missing characteristic ring system of beta-lactam antibiotics"

    # Look for amide groups (common in side chains)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide groups found in side chains"

    # Count atoms to ensure molecule is large enough to be an antibiotic
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 15:
        return False, "Molecule too small to be a beta-lactam antibiotic"

    # Count number of rings
    num_rings = len(Chem.GetSymmSSSR(mol))
    if num_rings < 2:
        return False, "Too few rings for a beta-lactam antibiotic"

    return True, "Contains beta-lactam ring and characteristic features of beta-lactam antibiotics"