"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Classifies a substance as a mineral based on its SMILES string using heuristic features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if classified as a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Common mineral elements
    mineral_elements = {'Ca', 'Mg', 'Fe', 'K', 'Na', 'Al', 'Si', 'S', 'Cl'}
    present_elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    
    # Check if key mineral elements are present
    if not any(element in present_elements for element in mineral_elements):
        return False, "No key mineral-forming elements found"
    
    # Common mineral-like ions
    ions = [
        Chem.MolFromSmarts('O=S(=O)([O-])[O-]'),  # sulfate
        Chem.MolFromSmarts('P(=O)([O-])([O-])[O-]'),  # phosphate
        Chem.MolFromSmarts('O=C([O-])[O-]'),  # carbonate
        Chem.MolFromSmarts('Cl'),  # chloride
    ]
    
    ion_presence = any(mol.HasSubstructMatch(ion) for ion in ions)
    if not ion_presence:
        return False, "No common mineral anion groups found"
    
    # Charge balance check (simple heuristic)
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, f"Charge imbalance detected: {total_charge}"

    # Check for crystal water (hydrates) pattern with multiple 'O'
    hydrates_pattern = Chem.MolFromSmarts('O.O.O')  # generalized for multiple water presence
    if mol.HasSubstructMatch(hydrates_pattern):
        return True, "Contains hydrate pattern typical in minerals"

    return True, "Contains mineral-forming elements and ions with charge balance"