"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:26811 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    Phosphatidyl-L-serine is an aminophospholipid with a phosphatidyl group esterified
    to the hydroxy group of L-serine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphatidyl-L-serine pattern
    ps_pattern = Chem.MolFromSmarts("[P](=[O])([O-])OC[C@@H](COC(=O)[C@@H](N)C(=O)O)OC(=O)[CX3]([CX4])[CX4]")
    if not mol.HasSubstructMatch(ps_pattern):
        return False, "Missing phosphatidyl-L-serine pattern"

    # Check for two fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    # Check molecular weight (typical range: 700-900 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700 or mol_wt > 900:
        return False, "Molecular weight outside typical range"

    return True, "Contains phosphatidyl group esterified to serine, with two fatty acid chains"