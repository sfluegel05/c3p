"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: secondary ammonium ion
Definition: An organic cation obtained by protonation of any secondary amino compound; 
           major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has an NHR2+ group where R represents organic groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            # Check formal charge
            if atom.GetFormalCharge() != 1:
                continue
                
            # Count attached hydrogens
            n_hydrogens = atom.GetTotalNumHs()
            
            # Get all neighbors
            neighbors = atom.GetNeighbors()
            
            # Count organic neighbors (carbon and other non-H atoms)
            organic_neighbors = []
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    organic_neighbors.append(neighbor)
                elif neighbor.GetAtomicNum() not in [1, 7, 8]:  # Not H, N, O
                    organic_neighbors.append(neighbor)
            
            # Secondary ammonium should have:
            # - Exactly 1 hydrogen
            # - Exactly 2 organic groups attached
            # - Formal charge of +1
            if n_hydrogens == 1 and len(organic_neighbors) == 2:
                # Verify the organic groups are not leaving groups (like OH)
                if all(len(list(n.GetNeighbors())) >= 2 for n in organic_neighbors):
                    return True, "Contains NHR2+ group (secondary ammonium ion)"
            
    return False, "No secondary ammonium ion found"

__metadata__ = {
    'chemical_class': {
        'name': 'secondary ammonium ion',
        'definition': 'An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.',
        'parents': ['organic cation', 'ammonium ion']
    }
}