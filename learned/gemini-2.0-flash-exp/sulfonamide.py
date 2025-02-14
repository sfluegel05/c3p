"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide has the structure R-S(=O)2-N-R', where R' is typically H or C.

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
    sulfonamide_core_pattern = Chem.MolFromSmarts("[S;X4](=[OX1])(=[OX1])[N;!H0]") #Sulfur is directly connected to 2 oxygens, the nitrogen must not have 0 H
    
    # Search for substructure match
    if not mol.HasSubstructMatch(sulfonamide_core_pattern):
         return False, "Sulfonamide core not found"
    
    #Verify that the nitrogen is attached to at least one C or H
    matches = mol.GetSubstructMatches(sulfonamide_core_pattern)
    for match in matches:
        nitrogen_atom_idx = match[2]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_atom_idx)
        n_neighbor_is_carbon_or_hydrogen = False
        for neighbor in nitrogen_atom.GetNeighbors():
            if neighbor.GetAtomicNum() in [1,6]:
                n_neighbor_is_carbon_or_hydrogen = True
                break
        if not n_neighbor_is_carbon_or_hydrogen:
            #Check if it is a sulfonylurea pattern (N-S(=O)2-N-C(=O)-N)
            sulfonylurea_pattern = Chem.MolFromSmarts("[N]-[S;X4](=[OX1])(=[OX1])-[N]-[CX3](=[OX1])-[N]")
            #Check if it is a sulfamide pattern N-S(=O)2-N
            sulfamide_pattern = Chem.MolFromSmarts("[N]-[S;X4](=[OX1])(=[OX1])-[N]")
            if not mol.HasSubstructMatch(sulfonylurea_pattern) and not mol.HasSubstructMatch(sulfamide_pattern):
                return False, "Nitrogen not bonded to at least one C or H and not part of a sulfonylurea/sulfamide"


    return True, "Contains a sulfonamide group (S(=O)2-N) with valid nitrogen"