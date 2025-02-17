"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol (A compound that contains exactly two free/alcoholic hydroxy groups)
Definition:
    “Diol” means that the molecule has exactly two free (alcoholic) –OH groups attached to carbon.
    –OH groups that are part of carbonyl-based functionalities (e.g. carboxylic acids) are ignored.
    
This classifier uses a SMARTS pattern [OX2H] to identify hydroxyl groups, then checks that each
matches the criteria.
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is classified as a diol (contains exactly two free/alcoholic hydroxy groups)
    based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES string and add explicit hydrogens to ensure –OH groups are present.
      2. Use a SMARTS pattern ([OX2H]) to find hydroxyl groups.
      3. For each –OH found, ensure that:
           - It is attached to at least one carbon.
           - The attached carbon is not double-bonded to another oxygen (i.e. it is not part of a carbonyl group).
      4. Count the valid –OH groups; if exactly two are found, the molecule is classified as a diol.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a diol (exactly two qualifying –OH groups), otherwise False.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to ensure hydroxyl hydrogens are present.
    mol = Chem.AddHs(mol)
    
    # Use SMARTS to identify hydroxyl groups [OX2H] (oxygen atom with two connections and one hydrogen)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    valid_oh_count = 0
    for match in oh_matches:
        oh_idx = match[0]
        oh_atom = mol.GetAtomWithIdx(oh_idx)
        # Check that the oxygen is attached to at least one carbon.
        neighbors = oh_atom.GetNeighbors()
        has_carbon = any(nbr.GetAtomicNum() == 6 for nbr in neighbors)
        if not has_carbon:
            continue  # skip if no carbon neighbor
            
        # For each neighboring carbon, ensure it is not part of a carbonyl (i.e. not double-bonded to another O)
        is_valid = False
        for nbr in neighbors:
            if nbr.GetAtomicNum() != 6:
                continue
            # Examine bonds around the carbon neighbor.
            carbon_atom = nbr
            has_carbonyl = False
            for bond in carbon_atom.GetBonds():
                other = bond.GetOtherAtom(carbon_atom)
                # Skip the bond with our hydroxyl oxygen.
                if other.GetIdx() == oh_idx:
                    continue
                # If any other bond is a double bond to an oxygen, mark this OH as non-alcoholic.
                if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_carbonyl = True
                    break
            if not has_carbonyl:
                is_valid = True
                break
        if is_valid:
            valid_oh_count += 1
            
    if valid_oh_count == 2:
        return True, "Molecule contains exactly two free (alcoholic) hydroxy groups and is classified as a diol."
    else:
        return False, f"Molecule contains {valid_oh_count} qualifying hydroxy groups, which does not match the diol definition (exactly two required)."

# Example usage:
# test_smiles = "OCCCCCCCCCCCCO"  # Example: docosane-1,3-diol
# result, reason = is_diol(test_smiles)
# print(result, reason)