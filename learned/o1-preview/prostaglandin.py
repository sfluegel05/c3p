"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is a naturally occurring compound derived from prostanoic acid,
    characterized by a C20 backbone with a cyclopentane ring and specific functional groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper function to count atoms of a particular atomic number
    def count_atoms(molecule, atomic_num):
        return sum(1 for atom in molecule.GetAtoms() if atom.GetAtomicNum() == atomic_num)
    
    # Check for cyclopentane ring
    cp_ring = Chem.MolFromSmarts("C1CCCC1")  # Simple cyclopentane ring
    if not mol.HasSubstructMatch(cp_ring):
        return False, "No cyclopentane ring found"
    
    # Check for carboxylic acid group (-C(=O)O)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
    
    # Count total number of carbons
    num_carbons = count_atoms(mol, 6)
    if num_carbons < 18 or num_carbons > 22:
        return False, f"Number of carbons is {num_carbons}, which is not in the expected range for prostaglandins (18-22)"
    
    # Check for hydroxyl groups attached to the cyclopentane ring
    # Prostaglandins often have hydroxyl groups on the cyclopentane ring
    cp_with_oh = Chem.MolFromSmarts("C1[C@H](O)C[C@H](O)C[C@H]1O")
    if mol.HasSubstructMatch(cp_with_oh):
        hydroxy_on_ring = True
    else:
        # Check for at least one hydroxyl group on the ring
        hydroxyl_on_ring = False
        cp_ring_matches = mol.GetSubstructMatches(cp_ring)
        for match in cp_ring_matches:
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0:
                        hydroxyl_on_ring = True
                        break
        if not hydroxyl_on_ring:
            return False, "No hydroxyl groups on cyclopentane ring found"
    
    # Check for side chains attached at appropriate positions
    # The side chains are usually aliphatic chains attached to the cyclopentane ring
    side_chain_pattern = Chem.MolFromSmarts("C1CCC(C1)-[C,C]")
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "Side chains attached to cyclopentane ring not found at expected positions"
    
    # Check for conjugated double bonds in side chains
    conjugated_double_bonds = Chem.MolFromSmarts("C=C-C=C")
    if not mol.HasSubstructMatch(conjugated_double_bonds):
        return False, "No conjugated double bonds in side chains"
    
    # If all checks pass, it is likely a prostaglandin
    return True, "Molecule matches prostaglandin structural features"

__metadata__ = {
    'chemical_class': {
        'name': 'prostaglandin',
        'definition': 'Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'error': '',
    'stdout': None
}