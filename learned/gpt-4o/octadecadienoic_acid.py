"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain fatty acid with exactly 18 carbons, two C=C double bonds,
    and a carboxylic acid group at one end.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool, str: True and reason if molecule is an octadecadienoic acid, False and reason otherwise.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use SMILES parsing to create a set of carbon chain to test for octadecadienoic acid
    # Find all chains as paths starting from a carboxylic acid group for proper context
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Missing carboxylic acid group"

    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    for match in carboxyl_matches:
        curr_atom = mol.GetAtomWithIdx(match[0])

        # Traverse a carbon chain starting from carboxylic carbon, ensuring 18 C and positional double bonds
        def dfs_find_carbon_chain(atom, depth, visited, max_depth, double_bond_count):
            visited[atom.GetIdx()] = True
            if atom.GetAtomicNum() == 6:
                if depth > max_depth[0]:
                    max_depth[0] = depth
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_bond_count += 1
                    next_atom = bond.GetOtherAtom(atom)
                    if not visited[next_atom.GetIdx()]:
                        dfs_find_carbon_chain(next_atom, depth + 1, visited, max_depth, double_bond_count)
            return double_bond_count
        
        visited = [False] * mol.GetNumAtoms()
        max_chain_length = [0]
        double_bond_count = dfs_find_carbon_chain(curr_atom, 0, visited, max_chain_length, 0)

        if max_chain_length[0] >= 17 and double_bond_count == 2:
            return True, "Contains an 18-carbon chain with 2 C=C double bonds and a carboxylic acid group"

    return False, "Does not meet structure requirements for an octadecadienoic acid"

# Sample usage
smiles_test = "CCCCCC\\C=C\\C=C/CCCCCCCC(O)=O"
result, reason = is_octadecadienoic_acid(smiles_test)
print(f"Result: {result}, Reason: {reason}")