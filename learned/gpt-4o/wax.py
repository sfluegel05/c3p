"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is considered a wax based on its SMILES string.
    A wax is characterized by long-chain hydrocarbons, typically containing esters,
    being malleable at ambient temperatures, favoring linear aliphatic chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a wax, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for ester linkage (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Finding and evaluating the longest carbon chains linked to ester oxygens
    carbon_chains = []
    
    # Find ester group atoms
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Explore each ester linkage to assess chain length
    for est_match in ester_matches:
        # Find longest carbon chain on both sides of the ester linkage
        c_chain_length = [0, 0]
        for i, atom_idx in enumerate(est_match[:-1]):  # Exclude the oxygen in COO
            explored = set()
            def explore_chain(a_idx, length=0):
                explored.add(a_idx)
                atom = mol.GetAtomWithIdx(a_idx)
                if atom.GetAtomicNum() == 6:  # Carbon
                    length += 1
                for neighbor in atom.GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if n_idx not in explored and neighbor.GetAtomicNum() == 6:
                        explore_chain(n_idx, length)
                return length
            
            c_chain_length[i] = explore_chain(atom_idx)
        
        carbon_chains.append(min(c_chain_length))  # Get the shorter length as both chains matter

    if len([chain for chain in carbon_chains if chain >= 14]) < 1:
        return False, "No sufficiently long carbon chains found for wax"
    
    # Calculate number of rotatable bonds to ascertain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, f"Insufficient rotatable bonds ({n_rotatable}), waxes should be flexible"
    
    # Check if the molecule contains rings, especially aromatic rings
    if mol.GetRingInfo().NumRings() > 0:
        # Still report as not a typical wax due to cyclic structures
        return False, "Contains non-linear structure components that are not typical for waxes"
    
    return True, f"Structure meets criteria for wax: Contains ester linkage, long carbon chains, and sufficient rotatable bonds ({n_rotatable})"