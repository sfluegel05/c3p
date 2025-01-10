"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:35277 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is an organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Get neighbors of the nitrogen atom
            neighbors = atom.GetNeighbors()
            carbon_count = 0
            other_neighbors = 0

            for neighbor in neighbors:
                atomic_num = neighbor.GetAtomicNum()
                if atomic_num == 6:
                    carbon_count += 1
                else:
                    other_neighbors += 1  # Nitrogen bonded to other heteroatom

            # Get total number of hydrogens (implicit + explicit)
            hydrogen_count = atom.GetTotalNumHs()
            
            # Nitrogen should have a total valence of 4 (tetrahedral)
            total_valence = carbon_count + other_neighbors + hydrogen_count

            # Check if nitrogen is bonded to exactly two carbons and hydrogens make up the rest
            if carbon_count == 2 and total_valence == 4:
                return True, "Contains secondary ammonium ion group (protonated secondary amine)"
        
    # If no matching nitrogen atom is found
    return False, "No secondary ammonium ion group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35277',
        'name': 'secondary ammonium ion',
        'definition': 'An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.',
        'parents': ['CHEBI:35274']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    }
}