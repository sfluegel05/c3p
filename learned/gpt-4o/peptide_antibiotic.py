"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 2:
        return False, "Insufficient peptide bonds found"

    # Check for cyclic peptide structures
    cycle_inds = Chem.GetSymmSSSR(mol)
    if not cycle_inds:
        return False, "No cyclic structures found, common in many peptide antibiotics"

    # Check sulfur presence, a common element in peptide antibiotics like daptomycin
    contains_sulfur = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())
    if contains_sulfur:
        return True, "Contains sulfur, which is common in peptide antibiotics"

    # Check for complex structures characteristic of peptide antibiotics
    heavy_atoms = mol.GetNumHeavyAtoms()
    num_rings = Chem.GetSSSR(mol)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    if heavy_atoms > 25 and num_rings > 1 and len(chiral_centers) > 3:
        return True, "Matches the complexity and structural features typical of peptide antibiotics"

    return False, "Structure does not match typical peptide antibiotics patterns"