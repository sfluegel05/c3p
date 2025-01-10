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
    A methyl sulfide is an aliphatic sulfide in which at least one of the groups attached to the sulfur is a methyl group.

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

    # Iterate over sulfur atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16 and atom.GetDegree() == 2:
            # Get neighbor atoms of sulfur
            neighbors = atom.GetNeighbors()
            # Check if both neighbors are carbons
            if all([nbr.GetAtomicNum() == 6 for nbr in neighbors]):
                # Check if at least one carbon is a methyl group
                for carbon_atom in neighbors:
                    # Exclude carbons in rings (aliphatic requirement)
                    if carbon_atom.IsInRing():
                        continue
                    # Get heavy atom neighbors of the carbon (exclude hydrogens and sulfur itself)
                    carbon_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != atom.GetIdx()]
                    if len(carbon_neighbors) == 0:
                        # Carbon is only connected to sulfur and hydrogens (methyl group)
                        return True, "Contains methyl sulfide group"
    return False, "No methyl sulfide group found"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXXXX',  # Replace with actual ChEBI ID if available
        'name': 'methyl sulfide',
        'definition': 'Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.',
        'parents': ['CHEBI:16683']  # Example parent class (aliphatic sulfide)
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        # Additional configuration parameters can be added here
    },
    # Additional metadata fields can be added as needed
}