"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:47926 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of exactly 3 amino acid residues connected by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bonds (-C(=O)-NH-)
    peptide_pattern = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])[#6]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    if len(peptide_matches) != 2:
        return False, f"Found {len(peptide_matches)} peptide bonds, need exactly 2"

    # Look for terminal amino group (-NH2)
    terminal_amine_pattern = Chem.MolFromSmarts("[NX3H2,NX3H3+][CX4H]")
    if not mol.HasSubstructMatch(terminal_amine_pattern):
        return False, "No terminal amino group found"

    # Look for terminal carboxyl group (-COOH or -COO-)
    terminal_carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(terminal_carboxyl_pattern):
        return False, "No terminal carboxyl group found"

    # Count nitrogens and oxygens to verify composition
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Tripeptides should have at least:
    # - 4 nitrogens (1 terminal NH2 + 2 peptide bonds + possible N in side chains)
    # - 5 oxygens (2 peptide C=O + 1 terminal COOH + possible O in side chains)
    if n_count < 3:
        return False, f"Too few nitrogens ({n_count}) for a tripeptide"
    if o_count < 5:
        return False, f"Too few oxygens ({o_count}) for a tripeptide"

    # Check molecular weight - should be in reasonable range for tripeptides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:  # Minimum weight for smallest possible tripeptide (e.g., Gly-Gly-Gly)
        return False, f"Molecular weight ({mol_wt}) too low for tripeptide"
    if mol_wt > 1000:  # Maximum weight for reasonable tripeptide
        return False, f"Molecular weight ({mol_wt}) too high for tripeptide"

    # Look for alpha carbon pattern typical in amino acids
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3H,NX3H2][CX4H]([#6])[CX3](=[OX1])")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if len(alpha_carbon_matches) < 2:  # Should find at least 2 alpha carbons
        return False, "Missing characteristic amino acid alpha carbon pattern"

    # If all checks pass, it's likely a tripeptide
    return True, "Contains 3 amino acid residues connected by 2 peptide bonds with terminal amino and carboxyl groups"