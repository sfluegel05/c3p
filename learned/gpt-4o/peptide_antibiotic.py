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

    # Look for peptide bond pattern (O=C-N)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "No peptide bonds found"
    
    # Check for cyclic peptide structures
    cycle_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[NX3][CX3](=O)@[NX3][CX3](=O)")
    if mol.HasSubstructMatch(cycle_pattern):
        return True, "Contains cyclic peptide bonds, a hallmark of many peptide antibiotics"

    # Check for sulfur present, common in peptide antibiotics like daptomycin
    contains_sulfur = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())
    if contains_sulfur:
        return True, "Contains sulfur, common in peptide antibiotics"

    # Assess overall complexity through various factors
    heavy_atoms = mol.GetNumHeavyAtoms()
    num_rings = Chem.GetSSSR(mol)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    # Adjusted threshold and complexity factors
    if heavy_atoms > 25 and num_rings >= 1 and len(chiral_centers) > 0:
        return True, "Has sufficient complexity in structure typical of peptide antibiotics"

    return False, "Structure does not match typical peptide antibiotics patterns"

# Note: Further specific pattern detection may improve the accuracy.