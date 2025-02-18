"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low for macromolecule: {mol_wt} Da"

    # Look for repeating units based on SMARTS patterns that might represent common monomeric units
    # Here, we use a simplistic approach by looking for repeating carbon chains or other simple units
    # This can be customized/expanded for specific types of macromolecules
    repeating_unit_pattern = Chem.MolFromSmarts("[*]~[*]~[*]~[*]~[*]")  # Very basic pattern
    if not mol.HasSubstructMatch(repeating_unit_pattern):
        return False, f"No repeating units found. Typical macromolecules have such units."

    return True, "Molecule is considered a macromolecule based on molecular weight and repeating units"