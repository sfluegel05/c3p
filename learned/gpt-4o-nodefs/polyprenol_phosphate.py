"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is characterized by a polyprenyl chain terminated by a phosphate or diphosphate group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Define SMARTS for polyprenyl chain with isoprene units
    isoprene_pattern = Chem.MolFromSmarts("C(=C)-C-C")
    
    # Find matching substructures for isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:  # Allow at least one isoprene unit to be more inclusive
        return False, "Too few isoprene units for polyprenol"

    # Check molecular weight typical for polyprenol phosphates
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for polyprenol phosphate"

    # Define SMARTS for phosphate and diphosphate groups (more focused)
    phosphate_pattern = Chem.MolFromSmarts("O=P(O)(O)[O-]")
    diphosphate_pattern = Chem.MolFromSmarts("O=P(O)(O)OP(=O)(O)[O-]")

    # Check if it contains a terminal phosphate or diphosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    
    if not (phosphate_matches or diphosphate_matches):
        return False, "No terminal phosphate or diphosphate group found"

    # Further validation for terminal group placement
    terminal_group_valid = any(mol.GetAtomWithIdx(match[-1]).GetDegree() == 1 for match in phosphate_matches + diphosphate_matches)
    if not terminal_group_valid:
        return False, "Terminal phosphate/diphosphate group not at the molecule end"

    return True, "Contains polyprenyl chain with terminal phosphate or diphosphate group"