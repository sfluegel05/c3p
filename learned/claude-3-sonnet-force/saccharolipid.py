"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is a lipid that contains a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify common carbohydrate moieties
    carbohydrate_patterns = [
        "[C@H]1[C@H]([C@@H]([C@H]([C@@H]1O)O)O)O",  # glucose
        "[C@H]1[C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O",  # galactose
        "[C@H]1[C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O",  # mannose
        "[C@H]1[C@@H]([C@H]([C@@H]([C@@H]1O)O)O)N(C)C=O",  # N-acetylglucosamine
        "[C@H]1[C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O[C@H]2[C@H]([C@@H]([C@H]([C@@H]2O)O)O)O"  # heptose
    ]
    carbohydrate_found = False
    for pattern in carbohydrate_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            carbohydrate_found = True
            break

    if not carbohydrate_found:
        return False, "No carbohydrate moiety found"

    # Identify common lipid moieties
    lipid_patterns = [
        "CCCCCCCCCCCC",  # fatty acid chain
        "C=CC(C)(C)C",  # isoprenoid
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)C(CCCCCCCCCCCCCC)C(=O)",  # mycolic acid
        "[C@@H]1[C@@H](O[C@H](CO)[C@@H](O)[C@@H]1OC(=O)[C@H](CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)"  # lipid A
    ]
    lipid_found = False
    for pattern in lipid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            lipid_found = True
            break

    if not lipid_found:
        return False, "No lipid moiety found"

    # Check for presence of both carbohydrate and lipid components
    carbohydrate_mcs = rdFMCS.FindMCS([mol, Chem.MolFromSmarts(carbohydrate_patterns[0])], matchValences=True, completeRingsOnly=True)
    lipid_mcs = rdFMCS.FindMCS([mol, Chem.MolFromSmarts(lipid_patterns[0])], matchValences=True, completeRingsOnly=True)

    if not carbohydrate_mcs.numAtoms or not lipid_mcs.numAtoms:
        return False, "No carbohydrate-lipid connectivity found"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for saccharolipid"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for saccharolipid"

    return True, "Contains both carbohydrate and lipid moieties"