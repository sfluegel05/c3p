"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:35457 beta-lactam antibiotic
An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.

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

    # Look for beta-lactam ring pattern
    beta_lactam_pattern = Chem.MolFromSmarts("[C&R1]1=[C&R2](N[C&R3]1[C&R4]=[O])")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Look for nitrogen-containing heterocycle
    heterocycle_pattern = Chem.MolFromSmarts("[&!r5&!r6]12[&!r5&!r6][&!r5&!r6][&!r5&!r6][&!r5&!r6]1[&!r5&!r6][&!r5&!r6][&!r5&!r6]2")
    heterocycle_match = mol.GetSubstructMatches(heterocycle_pattern)
    if not heterocycle_match:
        return False, "No nitrogen-containing heterocycle found"

    # Check for carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("[C&R5](=O)[O&R6]")
    carboxyl_match = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_match:
        return False, "No carboxyl group found"

    # Check for amine group
    amine_pattern = Chem.MolFromSmarts("[N&R7;!H0]")
    amine_match = mol.GetSubstructMatches(amine_pattern)
    if not amine_match:
        return False, "No amine group found"

    # If all conditions are met, classify as beta-lactam antibiotic
    return True, "Contains beta-lactam ring, nitrogen-containing heterocycle, carboxyl group, and amine group"