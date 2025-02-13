"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: CHEBI:24400 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a molecule with a carbohydrate (glycosyl) part linked to a lipid part via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycosidic bond pattern (C-O-C between carbohydrate and lipid)
    glycosidic_pattern = Chem.MolFromSmarts("[C,c][OX2][C,c]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond found"

    # Look for carbohydrate part (mono-, di-, or tri-saccharide)
    carbohydrate_pattern = Chem.MolFromSmarts("[C,c][OX2][C,c][OX2][C,c]")
    carbohydrate_matches = mol.GetSubstructMatches(carbohydrate_pattern)
    if len(carbohydrate_matches) == 0:
        return False, "No carbohydrate part found"

    # Look for lipid part (long carbon chain or sphingosine derivative)
    lipid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if len(lipid_matches) == 0:
        return False, "No lipid part found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be lipid"

    # Check molecular weight - glycolipids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycolipid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for glycolipid"
    if o_count < 5:
        return False, "Too few oxygens for glycolipid"

    return True, "Contains carbohydrate part linked to lipid part via glycosidic bond"