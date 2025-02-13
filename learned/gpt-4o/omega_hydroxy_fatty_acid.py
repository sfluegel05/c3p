"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl at the omega position.

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

    # Find carboxyl group at either end of a carbon chain
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group (-COOH) found"

    # Find hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl groups found"

    # Identify and verify the longest carbon chain
    longest_chain = None
    longest_length = 0
    for chain in rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True):
        c_chain_pattern = Chem.MolFromSmarts("C")
        if chain.HasSubstructMatch(c_chain_pattern):
            current_length = chain.GetNumAtoms()
            if current_length > longest_length:
                longest_length = current_length
                longest_chain = chain

    if longest_chain is None:
        return False, "No valid carbon chain found"

    # Ensure carboxyl group is at one terminal end
    matches = longest_chain.GetSubstructMatches(carboxyl_pattern)
    if len(matches) != 1:
        return False, "Incorrect number of terminal carboxyl groups"

    carboxyl_atom_index = matches[0][0]

    # Check if any hydroxyl group is at the opposite end of the carboxyl group in the longest chain
    opposite_terminal_hydroxyl = False
    for hydroxyl_match in hydroxyl_matches:
        hydroxyl_index = hydroxyl_match[0]
        hydroxyl_atom = longest_chain.GetAtomWithIdx(hydroxyl_index)
        for bond in hydroxyl_atom.GetBonds():
            bonded_atom = bond.GetOtherAtom(hydroxyl_atom)
            if bonded_atom.GetAtomicNum() == 6:  # Carbon
                carbon_index = bonded_atom.GetIdx()
                # Ensure this is a terminal carbon
                if carbon_index != carboxyl_atom_index and bonded_atom.GetDegree() == 1:
                    opposite_terminal_hydroxyl = True
                    break
    
    if not opposite_terminal_hydroxyl:
        return False, "No hydroxyl group at the omega position of the longest chain"

    return True, "Contains a terminal carboxyl group and an omega-position hydroxyl group, fitting the omega-hydroxy fatty acid structure"