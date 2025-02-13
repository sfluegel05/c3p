"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids with specific functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Simplified pattern focusing on a tetracyclic triterpenoid core
    tetracyclic_pattern = Chem.MolFromSmarts("C1CCC2C3C(CCC4C3CCC2C1)C4")
    if not mol.HasSubstructMatch(tetracyclic_pattern):
        return False, "No tetracyclic triterpenoid backbone found"

    # Check for at least two keto groups (C=O)
    keto_pattern = Chem.MolFromSmarts("[CX3]=O")
    keto_matches = mol.GetSubstructMatches(keto_pattern)
    if len(keto_matches) < 2:
        return False, f"Insufficient keto groups, found {len(keto_matches)}"

    # Check for multiple hydroxyl groups (-OH)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'H' for neighbor in atom.GetNeighbors()))
    if hydroxyl_count < 2:
        return False, f"Insufficient hydroxyl groups, found {hydroxyl_count}"

    # Cucurbitacins are typically larger molecules; use descriptive, not strict mass range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a typical cucurbitacin"

    return True, "Contains tetracyclic backbone with characteristic functional groups"