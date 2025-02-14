"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (-COOH) pattern
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group (-COOH) found"

    # Check for hydroxyl group (-OH) pattern
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group (-OH) found"

    # Check if hydroxyl is at the omega position
    num_atoms = mol.GetNumAtoms()
    oh_at_omega = False
    for match in hydroxyl_matches:
        oh_index = match[0]
        # Ensure the OH is not directly adjacent to the carboxyl group's C
        if mol.GetNumBonds() > 1:
            for bond in mol.GetBondWithIdx(oh_index).GetBeginAtom().GetBonds():
                neighbor = bond.GetOtherAtomIdx(oh_index)
                if neighbor == (num_atoms - 2):  # Avoid direct carboxylic C-OH adjacency
                    oh_at_omega = True
                    
    if not oh_at_omega:
        return False, "Hydroxyl group is not at the omega position"

    # Ensure the carbon chain is relatively linear and long
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 6:  # Typically omega-hydroxy fatty acids are longer
        return False, "Carbon chain is too short"

    return True, "Contains a terminal hydroxyl group and a carboxyl group, fitting the omega-hydroxy fatty acid structure"