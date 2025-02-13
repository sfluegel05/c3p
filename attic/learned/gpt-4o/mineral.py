"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Classifies a substance as a mineral based on its SMILES string using enhanced heuristic features.

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
    
    # Expanded list of common mineral elements
    mineral_elements = {
        'Ca', 'Mg', 'Fe', 'K', 'Na', 'Al', 'Si', 'S', 'Cl', 'O', 'Ni', 'Cu', 'Ba',
        'Zn', 'P', 'Cs', 'Sb', 'La', 'Pd'
    }
    present_elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    
    # Check if key mineral elements are present
    if not any(element in present_elements for element in mineral_elements):
        return False, "No key mineral-forming elements found"
    
    # Extended mineral-like ions including hydrate, silicate, and additional common anions
    ions = [
        Chem.MolFromSmarts('O=S(=O)([O-])[O-]'),  # sulfate
        Chem.MolFromSmarts('P(=O)([O-])([O-])[O-]'),  # phosphate
        Chem.MolFromSmarts('O=C([O-])[O-]'),  # carbonate
        Chem.MolFromSmarts('Cl'),  # chloride
        Chem.MolFromSmarts('F'),  # fluoride
        Chem.MolFromSmarts('[OH-]'),  # hydroxide
        Chem.MolFromSmarts('[Si](=O)([O-])[O-]')  # silicate
    ]
    
    ion_presence = any(mol.HasSubstructMatch(ion) for ion in ions)
    if not ion_presence:
        return False, "No common mineral anion groups found"
    
    # Charge balance check considering neutrality and possible net charges typical in minerals
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if not (total_charge == 0 or total_charge in {-1, 1, -2, 2}):  # Allow balanced and common mineral net charges
        return False, f"Charge inconsistency for recognizable mineral: {total_charge}"
    
    # Check for crystal water (hydrates) pattern with multiple 'O's
    hydrates_pattern = Chem.MolFromSmarts('O~O~O')  # look for a trio of connected oxygens as proxy for hydrate
    if mol.HasSubstructMatch(hydrates_pattern):
        return True, "Contains hydrate pattern typical in minerals"

    return True, "Contains mineral-forming elements and ions with charge balance"