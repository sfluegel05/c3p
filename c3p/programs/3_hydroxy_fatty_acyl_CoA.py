"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:28484 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is a molecule where the thiol group of coenzyme A is 
    condensed with the carboxyl group of a 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester bond (CoA linkage)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond (CoA linkage) found"

    # Check for 3-hydroxy group (hydroxy group on the third carbon of the fatty acid chain)
    # The pattern ensures the hydroxy group is on the third carbon from the thioester bond
    hydroxy_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4]([OX2H])[CX4]C(=O)[SX2]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3-hydroxy group found on the third carbon from the thioester bond"

    # Check for fatty acid chain (long carbon chain)
    # The pattern ensures the chain is connected to the thioester bond
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 3:
        return False, "Fatty acid chain too short or missing"

    # Check molecular weight - 3-hydroxy fatty acyl-CoA typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for 3-hydroxy fatty acyl-CoA"

    # Additional check to ensure the molecule is a 3-hydroxy fatty acyl-CoA
    # The pattern ensures the 3-hydroxy group is on the third carbon of the fatty acid chain
    specific_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4]([OX2H])[CX4]C(=O)[SX2]")
    if not mol.HasSubstructMatch(specific_pattern):
        return False, "Molecule does not match the specific pattern for 3-hydroxy fatty acyl-CoA"

    return True, "Contains thioester bond, 3-hydroxy group on the third carbon, and fatty acid chain"