"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    A semisynthetic derivative is defined as any organic molecular entity derived from a natural product by partial chemical synthesis.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for specific synthetic modifications
    acetylation_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")  # Acetylation
    halogenation_pattern = Chem.MolFromSmarts("[Cl,Br,I][CX4]")          # Halogenation
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")             # Ester group
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")             # Amide group
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")                     # Ether group
    sulfonamide_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])[NX3]") # Sulfonamide
    nitro_pattern = Chem.MolFromSmarts("[NX3](=[OX1])=[OX1]")            # Nitro group

    # Count the number of specific synthetic modifications
    synthetic_modifications = 0
    if mol.HasSubstructMatch(acetylation_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(halogenation_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(ester_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(amide_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(ether_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(sulfonamide_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(nitro_pattern):
        synthetic_modifications += 1

    # Check molecular weight - natural products are typically larger, but some derivatives can be smaller
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a natural product derivative"

    # Check for complex structure (e.g., multiple rings, chiral centers)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))

    # Check for natural product-like complexity
    is_complex = n_rings >= 2 and n_chiral_centers >= 1

    # If there are specific synthetic modifications and a complex structure, it is likely a semisynthetic derivative
    if synthetic_modifications >= 1 and is_complex:
        return True, "Contains specific synthetic modifications and a complex structure"

    # If no specific synthetic modifications are found, it is less likely to be a semisynthetic derivative
    return False, "No significant synthetic modifications or complex structure found"