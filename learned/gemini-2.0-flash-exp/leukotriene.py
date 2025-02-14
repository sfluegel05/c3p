"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is a C20 polyunsaturated fatty acid derivative with four double bonds,
    three of which are conjugated, derived from arachidonic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a leukotriene, False otherwise
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

   # 1. Check for a 20-carbon backbone with 4 double bonds and 3 conjugated.
    # This pattern ensures a specific backbone and the three conjugated double bond core.
    # It also looks for the carboxyl group to be part of the 20-carbon chain.
    # It requires the conjugated system to be within the 20 carbons.
    backbone_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]=[CX3,CX2]~[CX4,CX3]=[CX3,CX2]~[CX4,CX3]=[CX3,CX2]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3](=O)O")
    if not mol.HasSubstructMatch(backbone_pattern):
      return False, "Does not have a C20 backbone with the necessary conjugated double bond pattern"

    # 2. Check for four double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 4:
        return False, f"Less than 4 double bonds: found {len(double_bond_matches)}; not a leukotriene"

    # 3. Check for at least one hydroxyl
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"


    return True, "Meets the criteria for a leukotriene: 20 carbon backbone with four double bonds of which three are conjugated, and one hydroxyl and carboxylic acid."