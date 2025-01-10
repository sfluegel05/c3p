"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has an ester functional group with a decanoic acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Check for each ester if it can be classified as a decanoate ester
    for match in ester_matches:
        # Get the indices of the ester functional group
        carbonyl_c, ester_o = match[0], match[2]

        # Explore chains emerging from oxygen to find decanoic carbon chain
        # It should ideally be linear with 10 carbons starting from the ester oxygen
        atom = mol.GetAtomWithIdx(ester_o)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Verify it's a carbon
                carbon_chain = [neighbor.GetIdx()]
                current_atom = neighbor
                # Traverse the chain and count carbons
                while True:
                    for next_atom in current_atom.GetNeighbors():
                        if (next_atom.GetIdx() not in carbon_chain and 
                            next_atom.GetAtomicNum() == 6):
                            carbon_chain.append(next_atom.GetIdx())
                            current_atom = next_atom
                            break
                    else:
                        break

                if len(carbon_chain) == 10:
                    # We found a 10-carbon chain
                    return True, "Contains ester group with a decanoic acid (10-carbon) chain"

    return False, "Ester group without a proper decanoic acid chain"