"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are derived from prostanoic acid, a C20 fatty acid with a cyclopentane ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Core cyclopentane ring with at least one double bond
    core_pattern = Chem.MolFromSmarts("C1[C]([C])[C]([C])C[C]1=[C]")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Prostaglandin core ring not found"

    # 2. At least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"


    # 3. Check for two carbon chains attached to cyclopentane ring
    # simplified chain pattern checking for at least one chain of length 5+ carbons, including a carboxylic acid
    chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[C]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(chain_matches) < 2 or len(carboxylic_acid_matches) < 1:
       return False, "Not enough carbon chains or carboxylic acid group"

    # 4. Check for a minimum number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, "Too few carbons for a prostaglandin"

    # 5. Check for molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for prostaglandin"
    return True, "Meets prostaglandin criteria"