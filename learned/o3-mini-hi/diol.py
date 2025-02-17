"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol
Definition:
    A diol is defined as a compound that contains exactly two free (alcoholic) –OH groups 
    attached to sp3 carbons. Hydroxyl groups that are part of carbonyl-based functionalities 
    (e.g. in carboxylic acids) are ignored.
"""

from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is classified as a diol (contains exactly two free/alcoholic hydroxyl groups)
    based on its SMILES string.
    
    The function:
      - Converts the SMILES string into an RDKit molecule and adds explicit hydrogens.
      - Iterates over oxygen atoms to identify –OH groups that are:
            * Bound to exactly one hydrogen and one heavy atom.
            * The heavy atom is a carbon atom.
            * That carbon atom is sp3 hybridized.
            * The carbon is not involved in a double bond with another oxygen (which would indicate a carbonyl).
      - Counts these qualified –OH groups.
      - Returns True if exactly two such groups exist.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diol (exactly two qualifying –OH groups), otherwise False.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to ensure hydroxyl hydrogens are visible
    mol = Chem.AddHs(mol)
    
    valid_oh_indices = set()  # we will store indices of oxygen atoms that qualify as free OH groups
    
    # Iterate over all atoms looking for oxygen atoms that might be part of an alcoholic -OH
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:  # Only oxygen atoms
            continue

        # Check that the oxygen is involved in exactly 2 bonds
        # (it should be bound to one hydrogen and one heavy atom)
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            continue

        # There should be at least one hydrogen among its neighbors
        h_count = sum(1 for nb in neighbors if nb.GetAtomicNum() == 1)
        if h_count < 1:
            continue

        # Now, the other neighbor should be a carbon atom
        non_h_neighbors = [nb for nb in neighbors if nb.GetAtomicNum() != 1]
        if len(non_h_neighbors) != 1:
            continue  # skip if there is not exactly one heavy neighbor
        carbon = non_h_neighbors[0]
        if carbon.GetAtomicNum() != 6:
            continue  # we require that the -OH is attached to carbon

        # Check that the carbon is sp3 hybridized
        if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # Now ensure the carbon is not double-bonded to an oxygen (which would signal a carbonyl)
        is_carbonyl = False
        for bond in carbon.GetBonds():
            if bond.GetBondTypeAsDouble() == 2:  # is a double bond
                other = bond.GetOtherAtom(carbon)
                if other.GetAtomicNum() == 8:
                    is_carbonyl = True
                    break
        if is_carbonyl:
            continue

        # This oxygen qualifies as a free alcoholic hydroxyl group
        valid_oh_indices.add(atom.GetIdx())
    
    valid_oh_count = len(valid_oh_indices)
    
    if valid_oh_count == 2:
        return True, "Molecule contains exactly two free (alcoholic) hydroxyl groups attached to sp3 carbons and is classified as a diol."
    else:
        return False, f"Molecule contains {valid_oh_count} qualifying hydroxyl groups, which does not match the diol definition (exactly two required)."

# Example usage:
# test_smiles = "OCCCCCCCCCCCCO"  # Example: 1,12-dodecanediol
# result, reason = is_diol(test_smiles)
# print(result, reason)