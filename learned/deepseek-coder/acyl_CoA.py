"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:15347 acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a thioester bond (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond found"

    # Check for the presence of the CoA moiety
    # The CoA moiety includes a specific pattern of atoms and bonds, including the ADP group, pantothenic acid derivative, and cysteamine group.
    # This pattern captures the key features of CoA more accurately.
    coa_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][CX4][NX3][CX3](=[OX1])[CX4][CX4][SX2][CX3](=[OX1])[CX4]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for the presence of a carboxylic acid derivative attached to the thioester bond
    # This pattern captures the key features of a carboxylic acid derivative.
    carboxylic_acid_derivative_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0]")
    if not mol.HasSubstructMatch(carboxylic_acid_derivative_pattern):
        return False, "No carboxylic acid derivative found"

    # Ensure the thioester bond is connected to the CoA moiety
    thioester_match = mol.GetSubstructMatch(thioester_pattern)
    coa_match = mol.GetSubstructMatch(coa_pattern)
    if not any(atom in coa_match for atom in thioester_match):
        return False, "Thioester bond not connected to CoA moiety"

    # If all checks pass, the molecule is likely an acyl-CoA
    return True, "Contains a thioester bond, CoA moiety, and carboxylic acid derivative"

# Example usage:
# print(is_acyl_CoA("CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"))