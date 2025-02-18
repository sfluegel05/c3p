"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: CHEBI: ??? lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide consists of a peptide (multiple amide bonds) with an attached lipid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for peptide component: at least two amide bonds
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 2:
        return False, f"Found {len(amide_matches)} amide groups, need at least 2 for a peptide"

    # Check for lipid component: long chain (>=8 carbons) attached via amide/ester
    # Amide-linked lipid pattern: amide connected to at least 8 carbons
    amide_lipid_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])-!@*[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    # Ester-linked lipid pattern: ester connected to at least 8 carbons
    ester_lipid_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])-!@*[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")

    has_amide_lipid = mol.HasSubstructMatch(amide_lipid_pattern)
    has_ester_lipid = mol.HasSubstructMatch(ester_lipid_pattern)

    if not (has_amide_lipid or has_ester_lipid):
        return False, "No lipid chain (>=8 carbons) attached via amide or ester"

    # Optional: Verify lipid chain length using molecular weight or carbon count
    # Count total carbons in lipid chain (approximate)
    # This part is complex; relying on substructure match for simplicity

    return True, "Contains peptide with lipid chain attached via amide/ester"