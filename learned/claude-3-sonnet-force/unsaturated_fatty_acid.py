"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: CHEBI:36976 unsaturated fatty acid
Any fatty acid containing at least one C=C or C#C bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for unsaturation
    num_unsaturated_bonds = rdMolDescriptors.CalcNumUnsaturatedBonds(mol)
    if num_unsaturated_bonds == 0:
        return False, "No unsaturation found"
    
    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for carbon chain length
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 8 or n_carbons > 30:
        return False, "Carbon chain length outside typical range for fatty acids"
    
    # Check for disqualifying substructures
    disqualifying_patterns = [Chem.MolFromSmarts('c'), Chem.MolFromSmarts('N'), Chem.MolFromSmarts('O=C-O')]
    for pattern in disqualifying_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Molecule contains disqualifying substructure(s)"
    
    # Construct reason string
    reason = f"Contains {num_unsaturated_bonds} unsaturated bond(s), carboxylic acid group, and appropriate carbon chain length"
    
    return True, reason