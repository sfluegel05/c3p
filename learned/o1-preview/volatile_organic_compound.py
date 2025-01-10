"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: volatile organic compound
"""
from rdkit import Chem

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is any organic compound having an initial boiling point less than or equal to 250°C
    measured at a standard atmospheric pressure of 101.3 kPa.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    # Import necessary modules from RDKit
    from rdkit.Chem import Descriptors

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains carbon atoms (organic compound)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not an organic compound (contains no carbon atoms)"

    # Compute molecular descriptors
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1,6))
    num_rings = Chem.GetSSSR(mol)
    mol_weight = Descriptors.ExactMolWt(mol)

    # Identify key functional groups
    has_alcohol = mol.HasSubstructMatch(Chem.MolFromSmarts('[CX4][OX2H]'))  # primary or secondary alcohol
    has_carbonyl = mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3]=[OX1]'))  # carbonyl group
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)O[CX4]'))  # ester
    has_ether = mol.HasSubstructMatch(Chem.MolFromSmarts('[CX4][OX2][CX4]'))  # ether
    has_halogen = mol.HasSubstructMatch(Chem.MolFromSmarts('[F,Cl,Br,I]'))  # halide
    has_nitro = mol.HasSubstructMatch(Chem.MolFromSmarts('[$([NX3](=O)=O)]'))  # nitro group
    has_amine = mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]'))  # amine

    # Implement rule-based classification
    boiling_point = None
    reason = ""

    # Rules based on carbon count and functional groups
    if num_carbons <= 12:
        # Small molecules generally have boiling points ≤ 250°C
        boiling_point = "<= 250°C"
        reason = f"Contains {num_carbons} carbons; small organic molecule"
    elif has_alcohol and num_carbons <= 20:
        # Long-chain alcohols up to 20 carbons may have boiling points ≤ 250°C
        boiling_point = "<= 250°C"
        reason = f"Alcohol with {num_carbons} carbons"
    elif (has_carbonyl or has_ester or has_ether or has_halogen or has_nitro or has_amine) and num_carbons <= 20:
        # Molecules with polar functional groups
        boiling_point = "<= 250°C"
        reason = f"Contains polar functional group and {num_carbons} carbons"
    elif num_carbons <= 17 and num_rings > 0:
        # Cyclic compounds with moderate carbon counts
        boiling_point = "<= 250°C"
        reason = f"Cyclic compound with {num_carbons} carbons"
    else:
        # Assume boiling point > 250°C
        boiling_point = "> 250°C"
        reason = f"High molecular weight or lacks features indicative of low boiling point"

    # Final classification
    if boiling_point == "<= 250°C":
        return True, f"{reason}; Classified as VOC"
    else:
        return False, f"{reason}; Not classified as VOC"