"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid (SCFA) based on its SMILES string.
    An SCFA is defined as an aliphatic monocarboxylic acid with a chain length of C6 or less.
    Acceptable substituents are limited to hydroxy (-OH) and keto (=O) groups.
    If any other non-hydrocarbon substituents are present, the compound is not regarded as an SCFA.
    
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
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
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
    
    # Check for unacceptable functional groups
    # Allowed: carboxylic acid, hydroxy (-OH), keto (=O) groups attached to carbon chain
    # Disallowed: amino groups, multiple carboxylic acids, ester groups, ethers, etc.
    # SMARTS patterns for unacceptable groups
    disallowed_patterns = [
        Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]"),  # primary or secondary amines
        Chem.MolFromSmarts("C(=O)O[C,c]"),  # esters
        Chem.MolFromSmarts("C(=O)[O,N,S][!H]"),  # carboxylate salts, amides
        Chem.MolFromSmarts("[SX3](=O)"),  # sulfoxides
        Chem.MolFromSmarts("[PX3](=O)"),  # phosphine oxides
        Chem.MolFromSmarts("C(=O)C(=O)"),  # anhydrides
        Chem.MolFromSmarts("[OX2][OX2]"),  # peroxides
        Chem.MolFromSmarts("C#N"),  # nitriles
        Chem.MolFromSmarts("[F,Cl,Br,I]"),  # halogens
    ]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups"

    # Count the number of hydroxy and keto groups (excluding the one in carboxylic acid)
    # Hydroxy group attached to carbon (exclude carboxylic OH)
    hydroxy_pattern = Chem.MolFromSmarts("[C;!$(C=O)][OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    # Keto group (carbonyl not part of carboxylic acid)
    keto_pattern = Chem.MolFromSmarts("[C;$([!C](=O)[!O]);!$(C(=O)O)]")
    keto_matches = mol.GetSubstructMatches(keto_pattern)
    
    # Remove overlaps with carboxylic acid group
    carboxylic_acid_atoms = set()
    for match in carboxylic_acid_matches:
        carboxylic_acid_atoms.update(match)
    hydroxy_matches = [match for match in hydroxy_matches if all(idx not in carboxylic_acid_atoms for idx in match)]
    keto_matches = [match for match in keto_matches if all(idx not in carboxylic_acid_atoms for idx in match)]
    
    # Check for any other oxygen atoms not part of carboxylic acid, hydroxy, or keto groups
    allowed_oxygen_atoms = set()
    for match in hydroxy_matches + keto_matches:
        allowed_oxygen_atoms.update(match)
    allowed_oxygen_atoms.update([atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetIdx() in carboxylic_acid_atoms])
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetIdx() not in allowed_oxygen_atoms:
            return False, "Contains oxygen atoms in disallowed functional groups"

    # Validate the carbon chain length (longest chain including the carboxylic carbon)
    # Generate all possible paths starting from carboxylic acid carbon
    carboxylic_carbon_idx = carboxylic_acid_matches[0][0]  # The first atom in carboxylic acid match is the carbon
    lengths = []
    def dfs(atom_idx, visited):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 0
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited and neighbor.GetAtomicNum() == 6:
                length = dfs(n_idx, visited.copy())
                if length > max_length:
                    max_length = length
        return 1 + max_length  # Include current atom

    chain_length = dfs(carboxylic_carbon_idx, set())
    if chain_length > 6:
        return False, f"Carbon chain length is {chain_length}, which is more than 6"
    
    return True, "Molecule is an aliphatic monocarboxylic acid with chain length of C6 or less"