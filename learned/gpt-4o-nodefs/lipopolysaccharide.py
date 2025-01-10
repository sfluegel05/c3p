"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    The identification is based on generic features common in lipopolysaccharides like
    lipid chains linked to polysaccharide components and potential phosphate groups.
    This is, however, a rudimentary test due to the complexity and diversity of LPS structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule shows potential LPS features, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced lipid chain detection - includes ester and amide linkages
    lipid_pattern1 = Chem.MolFromSmarts("C(=O)O[C;R0][C;R0][C;R0][C;R0]")  # Ester-linked lipids
    lipid_pattern2 = Chem.MolFromSmarts("C(=O)N[C;R0][C;R0][C;R0][C;R0]")  # Amide-linked lipids
    lipid_matches1 = mol.HasSubstructMatch(lipid_pattern1)
    lipid_matches2 = mol.HasSubstructMatch(lipid_pattern2)

    if not (lipid_matches1 or lipid_matches2):
        return False, "No characteristic lipid chain motifs detected"

    # Enhanced saccharide pattern - looking for rings with multiple OH groups
    saccharide_pattern = Chem.MolFromSmarts("C1O[C@@H]([O])C(O)C(O)C1")  # Generic polyhydroxylated ring
    saccharide_matches = mol.GetSubstructMatches(saccharide_pattern)

    if len(saccharide_matches) < 1:
        return False, "No polysaccharide-like structures detected"

    # Check for phosphate groups often present in LPS
    phosphate_pattern = Chem.MolFromSmarts("[O-]P(=O)(O)O")
    phosphate_matches = mol.HasSubstructMatch(phosphate_pattern)
    
    if not phosphate_matches:
        return False, "No phosphate groups identified, which are common in LPS"

    # Check molecular size - though do not strictly enforce to avoid false negatives
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight is too low for typical LPS"

    return True, "Contains features indicating potential LPS: lipid chains, polysaccharides, and phosphate groups"