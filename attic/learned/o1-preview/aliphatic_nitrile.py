"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for nitrile group attached to aliphatic carbon
    aliphatic_nitrile_pattern = Chem.MolFromSmarts("[C;!a]#N")
    nitrile_matches = mol.GetSubstructMatches(aliphatic_nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group attached to aliphatic carbon found"
    
    # Check for aromatic rings in the molecule
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings > 0:
        return False, f"Molecule contains {num_aromatic_rings} aromatic ring(s), not purely aliphatic"
    
    # Check for aromatic atoms in the molecule
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if aromatic_atoms:
        return False, "Molecule contains aromatic atoms, not purely aliphatic"
    
    return True, "Contains nitrile group attached to aliphatic carbon and molecule is aliphatic"

__metadata__ = {   'chemical_class': {   'name': 'aliphatic nitrile',
                              'definition': 'Any nitrile derived from an aliphatic compound.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8}}