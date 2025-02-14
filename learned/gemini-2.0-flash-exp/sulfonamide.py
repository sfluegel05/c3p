"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide has the structure R-S(=O)2-N-R'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for the sulfonamide core
    sulfonamide_core_pattern = Chem.MolFromSmarts("S(=O)(=O)N")

    # Search for substructure match
    matches = mol.GetSubstructMatches(sulfonamide_core_pattern)
    if not matches:
         return False, "Sulfonamide core not found"

    # Verify bond order of S=O
    for match in matches:
       sulfur_atom_idx = match[0]
       sulfur_atom = mol.GetAtomWithIdx(sulfur_atom_idx)
       oxygen_count = 0
       for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), sulfur_atom_idx)
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    return False, "Sulfonamide S=O bond not double"
                oxygen_count += 1
       if oxygen_count != 2:
           return False, "Must have exactly 2 S=O bonds"


    # Check for at least one C or H directly attached to the N
    for match in matches:
        nitrogen_atom_idx = match[1]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_atom_idx)
        n_neighbor_is_carbon_or_hydrogen = False
        for neighbor in nitrogen_atom.GetNeighbors():
           if neighbor.GetAtomicNum() in [1, 6]:
              n_neighbor_is_carbon_or_hydrogen = True
              break # Found at least one carbon or hydrogen
        if not n_neighbor_is_carbon_or_hydrogen:
            return False, "Nitrogen atom not bonded to at least one C or H"

    return True, "Contains a sulfonamide group (S(=O)2-N)"