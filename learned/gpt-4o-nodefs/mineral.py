"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Common characteristics of minerals include metal cations paired 
    with anions or complex ions, often hydrated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for metal elements typically found in minerals
    metals = ['Ca', 'Mg', 'Fe', 'Ni', 'Ba', 'Zn', 'Cu', 'K', 'Cs', 'Na', 'Al', 'Sb', 'La']
    metal_found = any(atom.GetSymbol() in metals for atom in mol.GetAtoms())
    if not metal_found:
        return False, "No typical metal elements found"
    
    # Look for common inorganic anions like sulfate, phosphate, carbonate, etc.
    anions = ['[O-]S([O-])(=O)=O', '[O-]P([O-])(=O)=O', '[O-]C(=O)[O-]', '[Cl-]', '[OH-]', '[F-]']
    anion_found = any(mol.HasSubstructMatch(Chem.MolFromSmarts(anion)) for anion in anions)
    if not anion_found:
        return False, "No typical inorganic anions found"
    
    # Check for presence of water molecules indicating hydrates
    hydrate_pattern = Chem.MolFromSmarts("O")
    hydrate_count = len(mol.GetSubstructMatches(hydrate_pattern))
    
    return True, f"Contains metal elements and inorganic anions{' with hydration' if hydrate_count > 0 else ''}"