"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments, including flavones and flavonols,
    often with hydroxyl and glycosyl substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define generalized flavone core (2-phenylchromen-4-one) allowing substitutions
    flavone_core_smarts = '[O]=C1C=CC(=C2C=CC=CC2=O1)'

    # Define generalized flavonol core (3-hydroxyflavone) allowing substitutions
    flavonol_core_smarts = '[O]=C1C=CC(=C2C=CC=CC2=O1)O'

    # Create mol objects from SMARTS
    flavone_core_mol = Chem.MolFromSmarts(flavone_core_smarts)
    flavonol_core_mol = Chem.MolFromSmarts(flavonol_core_smarts)

    # Check for substructure match with flavone core
    if mol.HasSubstructMatch(flavone_core_mol):
        return True, "Molecule contains flavone core structure"

    # Check for substructure match with flavonol core
    if mol.HasSubstructMatch(flavonol_core_mol):
        return True, "Molecule contains flavonol core structure"

    # Check for O-glycosides of flavones/flavonols (attached sugars at hydroxyl positions)
    # Define SMARTS pattern for glycosylated flavonoids
    glycoside_pattern = Chem.MolFromSmarts('[O;!R][C@H]1O[C@H]([C@H]([C@H](O)[C@@H]1O)O)CO')

    # If molecule has a flavone core and glycoside
    if mol.HasSubstructMatch(flavone_core_mol) and mol.HasSubstructMatch(glycoside_pattern):
        return True, "Molecule is a glycosylated flavone (anthoxanthin)"

    # If molecule has a flavonol core and glycoside
    if mol.HasSubstructMatch(flavonol_core_mol) and mol.HasSubstructMatch(glycoside_pattern):
        return True, "Molecule is a glycosylated flavonol (anthoxanthin)"

    # Additional check for substitutions on flavone/flavonol cores
    # Allow for hydroxyl and methoxy substitutions on core
    substituents = Chem.MolFromSmarts('[OH,OC]')
    matches = mol.GetSubstructMatches(flavone_core_mol)
    if matches:
        # Check for substitutions attached to the core
        for match in matches:
            core_atoms = set(match)
            for atom in mol.GetAtoms():
                if atom.GetIdx() in core_atoms:
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() not in core_atoms:
                            if neighbor.GetSymbol() == 'O':
                                return True, "Molecule is a substituted flavone (anthoxanthin)"
    matches = mol.GetSubstructMatches(flavonol_core_mol)
    if matches:
        for match in matches:
            core_atoms = set(match)
            for atom in mol.GetAtoms():
                if atom.GetIdx() in core_atoms:
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() not in core_atoms:
                            if neighbor.GetSymbol() == 'O':
                                return True, "Molecule is a substituted flavonol (anthoxanthin)"

    # If no match is found
    return False, "Molecule does not match anthoxanthin criteria"