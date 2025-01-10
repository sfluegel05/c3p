"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine has a nitrogen atom bonded to two hydrogens and one carbon (NH2 group attached to a hydrocarbyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to accurately count hydrogen atoms
    mol = Chem.AddHs(mol)
    
    # Define primary amine pattern:
    # Nitrogen with valence 3, degree 1 (bonded to one carbon), two hydrogens,
    # excluding amides, imines, nitriles, nitro groups, and other nitrogen-containing functional groups
    primary_amine_smarts = "[NX3;H2;!$(N-[C,S]=[O,S,N]);!$(N#[C,N]);!$(N=*)]"
    primary_amine_pattern = Chem.MolFromSmarts(primary_amine_smarts)
    
    # Search for primary amine groups
    matches = mol.GetSubstructMatches(primary_amine_pattern)
    if not matches:
        return False, "No primary amine group found"
    
    for match in matches:
        nitrogen_idx = match[0]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Exclude nitrogens bonded to heteroatoms other than hydrogen or carbon
        bonded_atoms = [atom.GetAtomicNum() for atom in nitrogen_atom.GetNeighbors()]
        if any(atom_num not in [1, 6] for atom_num in bonded_atoms):
            continue  # Skip if bonded to atoms other than H or C
        
        return True, "Contains primary amine group (NH2 attached to hydrocarbyl group)"
    
    return False, "No valid primary amine group found"

__metadata__ = {
    'chemical_class': {
        'name': 'primary amine',
        'definition': 'A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.'
    },
    'config': {
        'additional_info': 'Classifies primary amines based on NH2 group attached to hydrocarbyl group.'
    }
}