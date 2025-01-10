"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Classifies a substance as a mineral based on its SMILES string using comprehensive analysis.
    
    Args:
        smiles (str): SMILES string of the substance
        
    Returns:
        bool: True if classified as a mineral, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated list of common elements found in minerals
    mineral_elements = {
        'Ca', 'Mg', 'Fe', 'K', 'Na', 'Al', 'Si', 'S', 'Cl', 'O', 'Ni', 'Cu', 'Ba', 
        'Zn', 'P', 'Cs', 'Sb', 'La', 'Pd', 'B'
    }
    present_elements = {atom.GetSymbol() for atom in mol.GetAtoms()}

    # Check if key mineral elements are present
    if not any(element in present_elements for element in mineral_elements):
        return False, "No key mineral-forming elements found"

    # Expanded patterns for mineral-like ions including borates
    ions = [
        Chem.MolFromSmarts('O=S(=O)([O-])[O-]'),  # sulfate
        Chem.MolFromSmarts('P(=O)([O-])([O-])[O-]'),  # phosphate
        Chem.MolFromSmarts('O=C([O-])[O-]'),  # carbonate
        Chem.MolFromSmarts('Cl'),  # chloride
        Chem.MolFromSmarts('F'),  # fluoride
        Chem.MolFromSmarts('[OH-]'),  # hydroxide
        Chem.MolFromSmarts('[Si](=O)([O-])[O-]'),  # silicate
        Chem.MolFromSmarts('B(O[H])O[H]'),  # simple borate
    ]

    ion_presence = any(mol.HasSubstructMatch(ion) for ion in ions)
    if not ion_presence:
        return False, "No common mineral anion groups found"
    
    # Charge balance check: Ensure neutrality or very simple net charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if not (total_charge == 0 or total_charge in {-1, 1}):
        return False, f"Charge inconsistency for recognizable mineral: {total_charge}"

    # Hydrates and other mineral characteristic patterns
    hydrates_pattern = Chem.MolFromSmarts('O~O~O')  # recognizing hydrate water pattern
    if mol.HasSubstructMatch(hydrates_pattern):
        return True, "Contains hydrate pattern typical in minerals"
    
    # If no specific hydrate pattern but essential characteristics match
    return True, "Contains mineral-forming elements and ions with charge balance"