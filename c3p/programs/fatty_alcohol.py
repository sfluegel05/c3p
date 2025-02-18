"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: CHEBI:30843 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a carbon chain length of 3 to 27 atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"

    # Ensure the molecule is aliphatic (no aromatic rings)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 0:
        return False, "Molecule contains aromatic rings, not aliphatic"

    # Count the number of carbon atoms in the longest aliphatic chain
    # Use the largest aliphatic chain to determine the carbon count
    aliphatic_chains = Chem.MolFromSmarts("[CX4H2]")
    matches = mol.GetSubstructMatches(aliphatic_chains)
    if not matches:
        return False, "No aliphatic chain found"
    
    # Find the longest aliphatic chain
    longest_chain = max(matches, key=len)
    c_count = len(longest_chain)
    if c_count < 3 or c_count > 27:
        return False, f"Carbon chain length {c_count} is outside the range of 3 to 27"

    # Count the number of oxygen atoms (should be at least one, but can be more)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found"

    # Check for unsaturation (double or triple bonds)
    unsaturation = sum(1 for bond in mol.GetBonds() if bond.GetBondType() in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE])
    if unsaturation > 0:
        # Allow unsaturated fatty alcohols
        pass

    # Check for other functional groups that are not typical in fatty alcohols
    # For example, carboxylic acids, esters, ketones, aldehydes, etc.
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")  # Carboxylic acid
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")  # Ester
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")  # Ketone
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")  # Aldehyde
    if (mol.HasSubstructMatch(carboxylic_acid_pattern) or 
        mol.HasSubstructMatch(ester_pattern) or 
        mol.HasSubstructMatch(ketone_pattern) or 
        mol.HasSubstructMatch(aldehyde_pattern)):
        return False, "Molecule contains non-aliphatic functional groups, not a fatty alcohol"

    # Ensure the molecule is primarily aliphatic
    # Check if the molecule has more than 3 non-aliphatic functional groups
    non_aliphatic_groups = Chem.MolFromSmarts("[!#6&!#1]")
    non_aliphatic_count = len(mol.GetSubstructMatches(non_aliphatic_groups))
    if non_aliphatic_count > 3:
        return False, "Molecule contains too many non-aliphatic functional groups"

    return True, f"Aliphatic alcohol with {c_count} carbon atoms and {o_count} oxygen atoms"