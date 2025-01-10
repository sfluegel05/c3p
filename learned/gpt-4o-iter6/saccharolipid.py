"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid should contain both a carbohydrate moiety and a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string provided."

    # Extended carbohydrate patterns including phosphorylated sugars, sulfonated sugars, and ketoses
    carb_patterns = [
        Chem.MolFromSmarts("OC1(CO)O[C@H](O)[C@H](O)[C@@H]1O"),                         # Simple sugar rings
        Chem.MolFromSmarts("OC1OC(CO)C(O)C(O)C1O"),                                   # Alternative sugar structures
        Chem.MolFromSmarts("O[C@@H]1[C@@H](O)[C@@H](O[C@H](COP(O)(O)=O)[C@@H]1O)O"),  # Phosphate sugar
        Chem.MolFromSmarts("O[C@@H]1O[C@H](O[C@H]2O[C@H](CO)[C@@H](O[C@H]2O)[C@H]2O[C@@H]([C@@H](O2)CO)O)[C@H]1O"),  # Sulfonated sugar
    ]

    # Broader lipid pattern: Long hydrocarbon chain with ester, amide or ketone linkage
    lipid_patterns = [
        Chem.MolFromSmarts("C(=O)OC[C@H1]CCCCCCCCCCCCCCCC"),    # Long chain esters
        Chem.MolFromSmarts("C(=O)N[C@H1]CCCCCCCCCCCCCCCC"),     # Long chain amides
        Chem.MolFromSmarts("C[C@H](CCCCCCCCCCCCCCC)C=O"),       # Ketone terminal carbonate
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)"),               # General long alkyl chain
    ]

    # Check for presence of a carbohydrate component
    carb_present = any(mol.HasSubstructMatch(carb_pattern) for carb_pattern in carb_patterns)
    if not carb_present:
        return False, "No adequate carbohydrate moiety found."

    # Check for presence of a lipid component
    lipid_present = any(mol.HasSubstructMatch(lipid_pattern) for lipid_pattern in lipid_patterns)
    if not lipid_present:
        return False, "No suitable long hydrocarbon chain (lipid component) found."

    return True, "Contains both carbohydrate and lipid components, indicating a saccharolipid."