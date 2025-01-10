"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Detect the presence of a sugar moiety
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H]1O[C@H]([C@H](O)[C@H](O)[C@H]1O)"),  # pyranose
        Chem.MolFromSmarts("[C@H]1O[C@H]([C@H](O)[C@H]1O)"),          # furanose
        Chem.MolFromSmarts("[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1"),  # alternative pyranose
    ]
    
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No sugar backbone found"
    
    # Look for sulfur replacements (S or -SR) specifically replacing oxygen
    sulfur_patterns = [
        Chem.MolFromSmarts("[C-S]"),  # Any carbon-sulfur bond, replaces hydroxyl
        Chem.MolFromSmarts("[O-S]"),  # Direct replacement of oxygen by sulfur
        Chem.MolFromSmarts("[S;D1]")  # Terminal sulfur suggesting replacement
    ]
    
    # Check that sulfur is correctly attached to sugar carbons
    correct_sulfur_attach = False
    for pattern in sulfur_patterns:
        matching_atoms = mol.GetSubstructMatches(pattern)
        for match in matching_atoms:
            # Verify sulfur attachment to sp3 carbon adjacent to sugar oxygens
            sulfur_atom = mol.GetAtomWithIdx(match[1])  # Getting the sulfur atom in [O-S]
            for neighbor in sulfur_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() >= 3:  # carbon directly linked
                    correct_sulfur_attach = True
                    break

    if not correct_sulfur_attach:
        return False, "No valid sulfur replacement detected"
    
    return True, "Contains carbohydrate structure with sulfur substitution(s)"