"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for methyl sulfide
    # Sulfur atom (not in ring), degree 2, bonded to a methyl group and an aliphatic carbon
    pattern = Chem.MolFromSmarts('[#16D2;!R]([CH3])[C;!R;!Ar]')
    
    # Find matches
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        for match in matches:
            sulfur_idx = match[0]
            sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
            # Ensure sulfur is only bonded to carbons (no heteroatoms)
            neighbor_atomic_nums = [nbr.GetAtomicNum() for nbr in sulfur_atom.GetNeighbors()]
            if all(num == 6 for num in neighbor_atomic_nums):
                return True, "Contains methyl sulfide group"
        return False, "Sulfur bonded to heteroatoms"
    else:
        return False, "No methyl sulfide group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16878',  # ChEBI ID for methyl sulfide
        'name': 'methyl sulfide',
        'definition': 'Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.',
        'parents': ['CHEBI:16683']  # Parent class (aliphatic sulfide)
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        # Additional configuration parameters can be added here
    },
    # Additional metadata fields can be added as needed
}