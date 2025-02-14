"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: CHEBI:32857 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid is a fatty acid with one or more alkyl substituents on the parent hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for alkyl branches (not just methyl groups)
    branch_pattern = Chem.MolFromSmarts("[CX4H2]")
    branch_matches = mol.GetSubstructMatches(branch_pattern)
    if len(branch_matches) < 1:
        return False, "No alkyl branches found"

    # Calculate the length of the longest carbon chain
    chain_length = AllChem.CalcChainLength(mol, AllChem.CalcPIFChain)
    if chain_length < 4:
        return False, "Carbon chain too short for fatty acid"

    # Check for allowed cyclic structures
    ring_info = mol.GetRingInfo()
    allowed_rings = [3, 4, 5, 6]  # Allow cyclopropane, cyclobutane, cyclopentane, cyclohexane
    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size not in allowed_rings:
            return False, f"Contains disallowed cyclic structure of size {ring_size}"

    # Check for halogen substituents
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    if mol.HasSubstructMatch(halogen_pattern):
        return False, "Contains halogen substituents"

    # Check for hydroxyl groups (common in branched-chain fatty acids)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H1]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        reason = "Contains a carboxylic acid group, alkyl branches, and hydroxyl group(s) on a long hydrocarbon chain"
    else:
        reason = "Contains a carboxylic acid group and alkyl branches on a long hydrocarbon chain"

    # Check molecular weight
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 150 or mol_weight > 500:
        return False, "Molecular weight outside typical range for branched-chain fatty acids"

    # All checks passed
    return True, reason