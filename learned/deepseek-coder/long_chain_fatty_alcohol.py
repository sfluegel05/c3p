"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol (C13 to C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol is defined as a fatty alcohol with a chain length ranging from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    # Find the longest carbon chain excluding functional groups
    def get_main_chain_length(mol):
        # Remove functional groups to focus on the carbon backbone
        simplified_mol = Chem.RWMol(mol)
        for atom in simplified_mol.GetAtoms():
            if atom.GetAtomicNum() != 6:  # Remove non-carbon atoms
                simplified_mol.RemoveAtom(atom.GetIdx())
        # Find the longest path in the simplified molecule
        return rdMolDescriptors.CalcLongestChain(simplified_mol)

    main_chain_length = get_main_chain_length(mol)
    if main_chain_length < 13 or main_chain_length > 22:
        return False, f"Main chain length is {main_chain_length}, must be between 13 and 22"

    # Check if any hydroxyl is attached to the main chain
    hydroxyl_on_main_chain = False
    for match in hydroxyl_matches:
        hydroxyl_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in hydroxyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                # Check if this carbon is part of the main chain
                if neighbor.GetDegree() >= 2:  # Main chain carbons typically have 2+ connections
                    hydroxyl_on_main_chain = True
                    break
        if hydroxyl_on_main_chain:
            break

    if not hydroxyl_on_main_chain:
        return False, "No hydroxyl group attached to the main carbon chain"

    # Allow some functional groups while still considering it a fatty alcohol
    allowed_patterns = [
        Chem.MolFromSmarts("[CX3]=[CX3]"),  # Double bonds
        Chem.MolFromSmarts("[CX2]#[CX2]"),  # Triple bonds
        Chem.MolFromSmarts("[OX2H]"),       # Hydroxyl groups
    ]

    # Check for truly disqualifying functional groups
    disallowed_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Carboxylic acid
        Chem.MolFromSmarts("[CX3](=O)[OX2H0]"),  # Ester
        Chem.MolFromSmarts("[NX3]"),             # Amines
        Chem.MolFromSmarts("[SX2]"),             # Sulfides
        Chem.MolFromSmarts("[CX3](=O)[NX3]"),    # Amides
    ]

    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups"

    return True, f"Contains a hydroxyl group attached to a carbon chain of length {main_chain_length}"