"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem

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
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group (-OH) found"

    # The omega hydroxyl should be furthest from the carboxyl group
    carboxyl_match = mol.GetSubstructMatch(carboxyl_pattern)  # Takes the first match
    carboxyl_carbon_index = carboxyl_match[0]

    # Gather indices of carbon atoms (excluding the one in the carboxyl group)
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_indices.remove(carboxyl_carbon_index)

    # Find the carbon atom bonded to the omega hydroxyl
    omega_hydroxyl_found = False
    for match in hydroxyl_matches:
        hydroxyl_oxygen_index = match[0]
        hydroxyl_atom = mol.GetAtomWithIdx(hydroxyl_oxygen_index)
        for bond in hydroxyl_atom.GetBonds():
            bonded_atom = bond.GetOtherAtom(hydroxyl_atom)
            if bonded_atom.GetAtomicNum() == 6:  # Check for bonding with carbon
                bonded_carbon_index = bonded_atom.GetIdx()
                
                # Check if this is the last carbon in the chain
                if bonded_carbon_index in carbon_indices and bonded_carbon_index != carboxyl_carbon_index:
                    omega_hydroxyl_found = True
                    break

    if not omega_hydroxyl_found:
        return False, "Hydroxyl group is not at the omega position"

    return True, "Contains a terminal hydroxyl group and a carboxyl group, fitting the omega-hydroxy fatty acid structure"