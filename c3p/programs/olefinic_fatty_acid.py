"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is characterized by a long hydrocarbon chain (usually ≥8 carbons)
    containing at least one carbon-carbon double bond (C=C), with a terminal carboxylic acid
    group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not typical for fatty acids"

    # Check for the presence of a terminal carboxylic acid group: [CX3](=O)[OX1H0-,OX2H1]
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No terminal carboxylic acid group found"
    
    # Check for the presence of at least one carbon-carbon double bond: [C]=[C]
    c_c_double_bond = Chem.MolFromSmarts("[C]=[C]")
    if not mol.HasSubstructMatch(c_c_double_bond):
        return False, "No carbon-carbon double bond found (no C=C)"

    # Count the total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, f"Too few carbon atoms ({c_count}), requires at least 8"

    # Optional: Ensure the carboxylic acid is at the end of a carbon chain
    # Find the carboxylic acid carbon
    carboxylic_atoms = mol.GetSubstructMatches(carboxylic_acid)
    carboxylic_carbons = [match[0] for match in carboxylic_atoms]
    # Check if carboxylic acid carbon is attached to a carbon chain
    is_terminal = False
    for idx in carboxylic_carbons:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(neighbors) == 1:
            is_terminal = True
            break
    if not is_terminal:
        return False, "Carboxylic acid group is not at the end of a carbon chain"

    return True, "Molecule is an olefinic fatty acid with at least one C=C double bond and sufficient carbon chain length"