"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: CHEBI:XXXXX epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a fatty acid containing an epoxide ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for epoxide ring (C1OC1)
    epoxide_pattern = Chem.MolFromSmarts("[C;r3][O;r3][C;r3]")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"

    # Check for a long carbon chain (at least 10 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10:
        return False, "Carbon chain too short to be a fatty acid"

    # Check molecular weight (epoxy fatty acids typically >200 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for epoxy fatty acid"

    return True, "Contains a fatty acid chain with an epoxide ring"