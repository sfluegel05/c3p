"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with an attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide backbone (at least 2 consecutive amide bonds)
    peptide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0][CX4][CX3](=[OX1])[NX3H0]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if len(peptide_matches) == 0:
        return False, "Not enough consecutive peptide bonds to form a peptide backbone"

    # Look for lipid moiety (long aliphatic chain with at least 8 carbons)
    lipid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if len(lipid_matches) == 0:
        return False, "No lipid moiety found (requires at least 8 carbon chain)"

    # Check if the lipid is attached to the peptide (amide or ester bond)
    lipid_attachment_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0,O][CX4]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    lipid_attachment_matches = mol.GetSubstructMatches(lipid_attachment_pattern)
    if len(lipid_attachment_matches) == 0:
        return False, "Lipid moiety not properly attached to peptide backbone"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be a lipid"

    # Check molecular weight - lipopeptides typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for lipopeptide"

    return True, "Contains peptide backbone with properly attached lipid moiety"