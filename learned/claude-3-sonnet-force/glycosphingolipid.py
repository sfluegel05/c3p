"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid that contains a carbohydrate residue attached to a ceramide backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ceramide backbone pattern
    ceramide_pattern = Chem.MolFromSmarts("[N;X3;H2,H1;!$(N(*)-*=[N,O,S])][CX3](=O)[CH0;X2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH0;X2]")
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No ceramide backbone found"

    # Check for carbohydrate residues
    carbohydrate_pattern = Chem.MolFromSmarts("O[C@H]1[C@H]([C@H]([C@@H]([C@H](O1)O)O)O)O")
    carbohydrate_matches = mol.GetSubstructMatches(carbohydrate_pattern)
    if not carbohydrate_matches:
        return False, "No carbohydrate residue found"

    # Check for glycosidic linkages
    glycosidic_pattern = Chem.MolFromSmarts("[O;X2][C@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic linkage found"

    # Check molecular weight and composition
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycosphingolipid"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 20 or o_count < 5:
        return False, "Insufficient carbon or oxygen content for glycosphingolipid"

    return True, "Contains a ceramide backbone with carbohydrate residue(s) attached via glycosidic linkage(s)"