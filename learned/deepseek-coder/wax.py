"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: Waxes
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is typically an ester of long-chain fatty acids and long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for long carbon chains attached to the ester group
    ester_atoms = ester_matches[0]
    ester_carbon = ester_atoms[1]  # The carbon in the ester group
    ester_oxygen = ester_atoms[0]  # The oxygen in the ester group

    # Get the atoms directly bonded to the ester carbon and oxygen
    ester_carbon_neighbors = [atom.GetIdx() for atom in mol.GetAtomWithIdx(ester_carbon).GetNeighbors()]
    ester_oxygen_neighbors = [atom.GetIdx() for atom in mol.GetAtomWithIdx(ester_oxygen).GetNeighbors()]

    # Remove the ester oxygen and carbon from the neighbors list
    ester_carbon_neighbors.remove(ester_oxygen)
    ester_oxygen_neighbors.remove(ester_carbon)

    # Check the length of the carbon chains attached to the ester group
    def get_chain_length(start_atom, exclude_atoms):
        visited = set()
        stack = [(start_atom, 0)]
        max_length = 0
        while stack:
            atom_idx, length = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            if length > max_length:
                max_length = length
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in exclude_atoms:
                    stack.append((neighbor.GetIdx(), length + 1))
        return max_length

    chain1_length = get_chain_length(ester_carbon_neighbors[0], {ester_carbon})
    chain2_length = get_chain_length(ester_oxygen_neighbors[0], {ester_oxygen})

    if chain1_length < 10 or chain2_length < 10:
        return False, "Chains too short to be considered a wax"

    # Check molecular weight - waxes typically have high molecular weights
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a wax"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for a wax"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (one ester group)"

    # Check for other functional groups that are not typical of waxes
    non_wax_functional_groups = ["[NX3]", "[SX2]", "[PX4]", "[CX3](=O)[OX2H1]", "[CX3](=O)[NX3]", "[CX3](=O)[OX2H0]"]
    for pattern in non_wax_functional_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Contains non-wax functional group: {pattern}"

    return True, "Contains a single ester group with long carbon chains, typical of waxes"