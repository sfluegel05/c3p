"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:46245 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    A ubiquinone is a benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone,
    typically with a polyprenoid side chain at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more accurate core structure pattern for 2,3-dimethoxy-5-methylbenzoquinone
    core_pattern = Chem.MolFromSmarts("[CH3]C1=C(C(=O)C(=C(C1=O)OC)OC)")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core 2,3-dimethoxy-5-methylbenzoquinone structure not found"

    # Check for the presence of a side chain at position 6
    # The side chain is typically a polyprenoid chain, which can vary in length
    # We look for a carbon chain attached to the core structure
    side_chain_pattern = Chem.MolFromSmarts("[CH3]C1=C(C(=O)C(=C(C1=O)OC)OC)[!H]")
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "No side chain found at position 6"

    # Count the number of rotatable bonds to ensure the presence of a flexible side chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Side chain too short or rigid to be a polyprenoid chain"

    # Check molecular weight - ubiquinones typically have a higher molecular weight due to the side chain
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for ubiquinone"

    # Count carbons and oxygens to ensure the presence of the core structure and side chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for ubiquinone"
    if o_count < 4:
        return False, "Too few oxygens for ubiquinone"

    return True, "Contains core 2,3-dimethoxy-5-methylbenzoquinone structure with a polyprenoid side chain at position 6"