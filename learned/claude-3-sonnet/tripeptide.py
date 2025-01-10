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

    # Look for peptide bonds (-C(=O)-N-) with more flexible pattern
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    # For a tripeptide, we need exactly 2 peptide bonds
    if len(peptide_matches) < 2:
        return False, f"Found {len(peptide_matches)} peptide bonds, need at least 2"
    if len(peptide_matches) > 3:
        return False, f"Found {len(peptide_matches)} peptide bonds, too many for a tripeptide"

    # Look for amino acid residue pattern - more flexible version
    aa_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])")
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    if len(aa_matches) < 2:
        return False, "Too few amino acid residue patterns found"

    # Check for terminal groups with more flexible patterns
    terminal_amine_pattern = Chem.MolFromSmarts("[NX3H2,NX3H3+,$(N([H])[C])]")
    terminal_carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-]")
    
    if not (mol.HasSubstructMatch(terminal_amine_pattern) and 
            mol.HasSubstructMatch(terminal_carboxyl_pattern)):
        return False, "Missing terminal amino or carboxyl group"

    # Count nitrogens and oxygens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Revised minimum requirements:
    # - At least 3 nitrogens (1 terminal + 2 peptide bonds)
    # - At least 4 oxygens (2 peptide C=O + 1 terminal COOH)
    if n_count < 3:
        return False, f"Too few nitrogens ({n_count}) for a tripeptide"
    if o_count < 4:
        return False, f"Too few oxygens ({o_count}) for a tripeptide"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180:  # Adjusted minimum weight
        return False, f"Molecular weight ({mol_wt}) too low for tripeptide"
    if mol_wt > 1200:  # Adjusted maximum weight
        return False, f"Molecular weight ({mol_wt}) too high for tripeptide"

    # Look for backbone connectivity
    backbone_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])[NX3][CX4][CX3](=[OX1])[NX3][CX4][CX3](=[OX1])")
    if not mol.HasSubstructMatch(backbone_pattern):
        # Try alternative pattern for cyclic residues (like proline)
        alt_backbone = Chem.MolFromSmarts("[$(N([C])[C])[CX4][CX3](=[OX1])[$(N([C])[C])][CX4][CX3](=[OX1])[$(N([C])[C])][CX4][CX3](=[OX1])")
        if not mol.HasSubstructMatch(alt_backbone):
            return False, "Missing characteristic tripeptide backbone pattern"

    return True, "Contains 3 amino acid residues connected by peptide bonds"