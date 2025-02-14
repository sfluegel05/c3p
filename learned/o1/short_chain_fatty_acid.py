"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid (SCFA) based on its SMILES string.
    An SCFA is defined as an aliphatic monocarboxylic acid with a chain length of 6 or less carbons.
    If any non-hydrocarbon substituent is present (other than allowed hydroxy or keto groups),
    the compound is not regarded as an SCFA.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"
    
    # Ensure molecule is aliphatic (no aromatic atoms)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic atoms"
    
    # Ensure all atoms are C, H, or O
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 8]:
            return False, f"Contains atom other than C, H, O: {atom.GetSymbol()}"
    
    # Check the total number of oxygen atoms
    num_oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygen_atoms > 3:
        return False, f"Contains {num_oxygen_atoms} oxygen atoms, maximum allowed is 3"
    
    # Check for allowed functional groups outside carboxylic acid
    # Allowed additional functional groups: one hydroxy (-OH) or one keto (=O)
    # Find hydroxyl groups (excluding carboxylic acid OH)
    hydroxy_pattern = Chem.MolFromSmarts("[C][O;H1]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    # Exclude the OH in the carboxylic acid
    carboxylic_acid_atoms = set()
    for match in carboxylic_acid_matches:
        carboxylic_acid_atoms.update(match)
    hydroxy_matches = [match for match in hydroxy_matches if match[1] not in carboxylic_acid_atoms]
    
    # Find keto groups (excluding the carbonyl in carboxylic acid)
    keto_pattern = Chem.MolFromSmarts("C(=O)[C!O]")
    keto_matches = mol.GetSubstructMatches(keto_pattern)
    keto_matches = [match for match in keto_matches if match[0] not in carboxylic_acid_atoms]
    
    total_additional_oxygen_atoms = len(hydroxy_matches) + len(keto_matches)
    if total_additional_oxygen_atoms > 1:
        return False, f"Contains {total_additional_oxygen_atoms} additional hydroxy/keto groups, maximum allowed is 1"
    
    # Check for other disallowed functional groups
    # Disallowed: amines, halogens, esters, ethers, nitriles, multiple hydroxy or keto groups, etc.
    disallowed_patterns = [
        Chem.MolFromSmarts("[NX3;!$(N-C=O)]"),  # amines excluding amides
        Chem.MolFromSmarts("[SX2]"),  # thiols
        Chem.MolFromSmarts("[SX3](=O)"),  # sulfoxides
        Chem.MolFromSmarts("[SX4](=O)(=O)"),  # sulfones
        Chem.MolFromSmarts("C(=O)O[C,c]"),  # esters
        Chem.MolFromSmarts("[O;!H0;!$([O][C]=O)]([C,c])[C,c]"),  # ethers
        Chem.MolFromSmarts("C#N"),  # nitriles
        Chem.MolFromSmarts("[F,Cl,Br,I]"),  # halogens
        Chem.MolFromSmarts("[P]"),  # phosphorous
        Chem.MolFromSmarts("[S]"),  # sulfur
    ]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups"
    
    # Validate the carbon chain length (longest chain including the carboxylic carbon)
    # Identify the carboxylic acid carbon index
    carboxylic_carbon_idx = carboxylic_acid_matches[0][0]  # The first atom in carboxylic acid match is the carbon
    
    # Find all paths starting from carboxylic carbon
    # We will use a depth-first search to find the longest carbon chain
    visited = set()
    max_chain_length = [0]  # use list for mutable integer

    def dfs(node_idx, length):
        visited.add(node_idx)
        atom = mol.GetAtomWithIdx(node_idx)
        if atom.GetAtomicNum() != 6:  # Only consider carbon atoms
            return
        if length > max_chain_length[0]:
            max_chain_length[0] = length
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(node_idx, n_idx)
            if n_idx not in visited and bond is not None:
                n_atom = mol.GetAtomWithIdx(n_idx)
                if n_atom.GetAtomicNum() == 6:
                    dfs(n_idx, length + 1)
        visited.remove(node_idx)

    dfs(carboxylic_carbon_idx, 1)

    if max_chain_length[0] > 6:
        return False, f"Carbon chain length is {max_chain_length[0]}, which is more than 6"
    
    return True, "Molecule is an aliphatic monocarboxylic acid with chain length of 6 or less"