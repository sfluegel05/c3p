from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyridinedicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a pyridinedicarboxylic acid.
    These are pyridines with exactly two carboxylic acid groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyridinedicarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for pyridine core
    pyridine_pattern = Chem.MolFromSmarts('n1ccccc1')
    if not mol.HasSubstructMatch(pyridine_pattern):
        return False, "No pyridine core found"

    # Check for carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)O')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)

    if len(carboxyl_matches) != 2:
        return False, f"Found {len(carboxyl_matches)} carboxylic acid groups, need exactly 2"

    # Verify carboxylic acids are directly attached to pyridine ring
    pyridine_atoms = mol.GetSubstructMatch(pyridine_pattern)
    if not pyridine_atoms:
        return False, "Could not map pyridine atoms"

    pyridine_carbons = set([atom_idx for atom_idx in pyridine_atoms 
                           if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C'])

    for match in carboxyl_matches:
        carboxyl_carbon = match[0]  # First atom in SMARTS pattern is the carbon
        # Get atom attached to carboxyl group
        attached_atoms = [x.GetIdx() for x in mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors() 
                         if x.GetIdx() not in match]
        if not any(atom_idx in pyridine_carbons for atom_idx in attached_atoms):
            return False, "Carboxylic acid group not directly attached to pyridine ring"

    # Get positions of carboxylic acids
    positions = []
    ring_atoms = list(pyridine_atoms)
    for match in carboxyl_matches:
        carboxyl_carbon = match[0]
        for neighbor in mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors():
            if neighbor.GetIdx() in ring_atoms:
                positions.append(ring_atoms.index(neighbor.GetIdx()) + 1)
    positions.sort()

    return True, f"Pyridinedicarboxylic acid with COOH groups at positions {','.join(map(str,positions))}"
# Pr=1.0
# Recall=1.0