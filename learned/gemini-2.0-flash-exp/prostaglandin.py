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

    # 1. Specific cyclopentane ring with substitutions
    # core structure with double bond and functional groups.
    core_pattern = Chem.MolFromSmarts("[C]1([C]([O])([H])[C]([O])([H])[C]([H])([H])[C]([H])([H])1)~[C]~[C](=[O])O")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Prostaglandin core pattern not found"

    # 2. Check for at least one double bond
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond present"

    # 3. Check for presence of 2 carbon chains attached to cyclopentane ring (at least one carboxyl)
    chain_pattern_1 = Chem.MolFromSmarts("[C]1([C]([O])[H])([C]([O])[H])[C]([H])[C]([H])1~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    chain_pattern_2 = Chem.MolFromSmarts("[C]1([C]([O])[H])([C]([O])[H])[C]([H])[C]([H])1~[C]~[C]~[C]~[C]~[C]")
    if not (mol.HasSubstructMatch(chain_pattern_1) or mol.HasSubstructMatch(chain_pattern_2)):
        return False, "Prostaglandin chains not found"

    
    # 4 Check for two hydroxyl groups attached to ring
    alcohol_pattern = Chem.MolFromSmarts("[C]1([C]([O])[H])([C]([O])[H])[C][C]1")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "Did not find 2 hydroxyl groups attached to the ring"
    
    #5 Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
         return False, "No carboxylic acid group"

    return True, "Meets prostaglandin criteria"