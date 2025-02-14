"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:15948 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is any nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of a ribose sugar
    ribose_pattern = Chem.MolFromSmarts("[C@H]1([C@H](CO)O[C@@H](O)[C@H]1O)O")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar found"

    # Check for presence of a nucleobase component
    # Using a general aromatic ring pattern to account for non-standard bases
    nuc_base_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(nuc_base_pattern):
        return False, "No nucleobase component found"

    # Check for any bond between the ribose and the nucleobase
    # This allows for different types of glycosidic bonds or linkages
    bonds = mol.GetSubstructMatches(ribose_pattern)[0]
    ribose_atoms = set([mol.GetBondWithIdx(bond).GetBeginAtomIdx() for bond in bonds])
    base_atoms = set([atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8, 17]])  # Exclude C, H, O, Cl
    sugar_base_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts("[C,N,O,P]~[C,N,O,P]"))
    has_sugar_base_bond = any(set([bond[0], bond[1]]).intersection(ribose_atoms) and set([bond[0], bond[1]]).intersection(base_atoms) for bond in sugar_base_bonds)

    if not has_sugar_base_bond:
        return False, "No bond found between ribose sugar and nucleobase"

    # Additional checks or heuristics can be added here

    return True, "Contains a ribose sugar and a nucleobase component with a bond between them"