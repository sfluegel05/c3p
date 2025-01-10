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
    
    # Exclude molecules containing metals or disallowed elements
    allowed_atomic_nums = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}  # H, C, N, O, F, P, S, Cl, Br, I
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atomic_nums:
            return False, f"Molecule contains disallowed element with atomic number {atomic_num}"
    
    # Find nitrile groups: carbon triple-bonded to nitrogen
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group found"
    
    # Check that nitrile carbon is attached to an aliphatic (non-aromatic) carbon
    for match in nitrile_matches:
        nitrile_c_idx = match[0]
        nitrile_n_idx = match[1]
        nitrile_c_atom = mol.GetAtomWithIdx(nitrile_c_idx)
        nitrile_n_atom = mol.GetAtomWithIdx(nitrile_n_idx)
        
        # Get neighbors of nitrile carbon, excluding the nitrogen
        neighbors = [nbr for nbr in nitrile_c_atom.GetNeighbors() if nbr.GetIdx() != nitrile_n_idx]
        
        if not neighbors:
            return False, "Nitrile carbon has no attached group"
        
        # Check that the nitrile carbon is attached to a non-aromatic carbon
        attached_atom = neighbors[0]
        if attached_atom.GetAtomicNum() != 6:
            return False, "Nitrile carbon is not bonded to a carbon atom"
        if attached_atom.GetIsAromatic():
            return False, "Nitrile carbon is attached to an aromatic carbon"
        
        # Optional: check that the attached carbon is aliphatic (sp3 hybridized)
        if attached_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            return False, "Nitrile carbon is not attached to an aliphatic carbon"
    
    # If passed all checks, it's an aliphatic nitrile
    return True, "Contains nitrile group attached to an aliphatic carbon"

__metadata__ = {   
    'chemical_class': {   
        'name': 'aliphatic nitrile',
        'definition': 'Any nitrile derived from an aliphatic compound.'
    },
    'config': {   
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8
    }
}