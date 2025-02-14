"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:28304 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is a saponin derived from a hydroxysteroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone patterns
    steroid_patterns = [
        Chem.MolFromSmarts("[C@H]1CC[C@]2(C)[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@@]21C"),  # Androstane
        Chem.MolFromSmarts("[C@H]1CC[C@]2(C)[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@@]12C"),  # Pregnane
        Chem.MolFromSmarts("[C@H]1CC[C@]2(C)[C@@H]3C[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@@]12C"),  # Cholane
        # Add more steroid backbone patterns as needed
    ]
    has_steroid_backbone = any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns)
    if not has_steroid_backbone:
        return False, "No steroid backbone found"
    
    # Look for glycoside groups (-O-C(=O)-C-C-O-)
    glycoside_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4][CX4][OX2]")
    glycoside_matches = mol.GetSubstructMatches(glycoside_pattern)
    if not glycoside_matches:
        return False, "No glycoside groups found"
    
    # Check for additional structural features
    has_other_rings = bool(Chem.GetSymmSSSR(mol))
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Apply additional rules or heuristics as needed
    if mol_wt < 500:
        return False, "Molecular weight too low for steroid saponin"
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for steroid saponin"
    if not has_other_rings:
        return False, "No additional ring systems found"
    
    return True, "Contains steroid backbone and glycoside groups, with additional structural features consistent with steroid saponins"