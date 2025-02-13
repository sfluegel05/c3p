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

    # Check for the presence of multiple peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 5:
        return False, "Insufficient peptide bonds found"

    # Check for unique amino acids or cyclic structures patterns
    cyclic_pattern = Chem.MolFromSmarts("C1CCC(CC1)C1=CNC2=NC=CC=C2C1")
    if mol.HasSubstructMatch(cyclic_pattern):
        return True, "Contains unique cyclic patterns often found in peptide antibiotics"

    # Check sulfur or specific thiazole/oxazole presence, indicative of some antibiotics
    thiazole_pattern = Chem.MolFromSmarts("c1sccc1")
    oxazole_pattern = Chem.MolFromSmarts("c1noc1")
    if mol.HasSubstructMatch(thiazole_pattern) or mol.HasSubstructMatch(oxazole_pattern):
        return True, "Contains thiazole/oxazole rings, which are common in some peptide antibiotics"

    # Check overall complexity: large molecule with multiple rings and chiral centers
    heavy_atoms = mol.GetNumHeavyAtoms()
    num_rings = Chem.GetSSSR(mol)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    if heavy_atoms > 25 and num_rings > 3 and len(chiral_centers) > 5:
        return True, "Complex structure fits the profile of many peptide antibiotics"

    return False, "Structure does not match typical peptide antibiotics patterns"