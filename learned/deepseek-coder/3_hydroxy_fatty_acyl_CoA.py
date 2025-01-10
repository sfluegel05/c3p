"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:28494 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is a CoA ester of a 3-hydroxy fatty acid.

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

    # Look for CoA moiety pattern (thiol group attached to adenine, ribose, and phosphate groups)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Look for 3-hydroxy fatty acid chain (hydroxyl at third position from carboxyl group)
    # Flexible pattern to account for double bonds and branching
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]([OH])~[CX4,CX3]~[CX4,CX3](=O)S")
    matches = mol.GetSubstructMatches(hydroxy_fatty_acid_pattern)
    if not matches:
        return False, "No 3-hydroxy fatty acid chain found"

    # Validate the position of the hydroxyl group relative to the carboxyl group
    # The hydroxyl group should be at the third position from the carboxyl group
    for match in matches:
        # Get the indices of the matched atoms
        hydroxyl_index = match[2]  # Index of the hydroxyl group
        carboxyl_index = match[4]  # Index of the carboxyl carbon
        # Check if the hydroxyl group is at the third position from the carboxyl group
        if mol.GetBondBetweenAtoms(hydroxyl_index, carboxyl_index) is not None:
            return True, "Contains CoA moiety with 3-hydroxy fatty acid chain attached via ester bond"

    return False, "Hydroxyl group not at the third position from the carboxyl group"