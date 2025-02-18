"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule can be classified as a semisynthetic derivative based on broad chemical patterns.
    Semisynthetic derivatives are typically derived from natural products with additional synthetic modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if evidence suggests a semisynthetic derivative, False if not a confident classification
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for presence of patterns commonly introduced in semisynthetic derivatives
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")  # Ester linkage
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[#6]")  # Amide linkage
    halogen_pattern = Chem.MolFromSmarts("[#6][F,Cl,Br,I]")  # Halogen bonded to carbon

    ester_matches = mol.GetSubstructMatches(ester_pattern)
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    halogen_matches = mol.GetSubstructMatches(halogen_pattern)

    # Determine if there are significant matches to provide indirect evidence
    if len(ester_matches) > 0 or len(amide_matches) > 0 or len(halogen_matches) > 0:
        reason = "Molecule contains synthetic-like modifications: "
        reason += f"Ester linkages ({len(ester_matches)}), " if len(ester_matches) > 0 else ""
        reason += f"Amide linkages ({len(amide_matches)}), " if len(amide_matches) > 0 else ""
        reason += f"Halogen substitutions ({len(halogen_matches)})" if len(halogen_matches) > 0 else ""
        return True, reason.strip(", ")
    
    return False, "No significant synthetic patterns detected, classification not confident"