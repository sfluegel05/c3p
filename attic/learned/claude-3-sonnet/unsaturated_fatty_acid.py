"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: unsaturated fatty acid
Definition: Any fatty acid containing at least one C=C or C#C bond
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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
    
    # Check for carboxylic acid group (-COOH)
    carboxylic_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for double bonds (C=C)
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    
    # Check for triple bonds (C#C)
    triple_bond_pattern = Chem.MolFromSmarts('C#C')
    triple_bonds = mol.GetSubstructMatches(triple_bond_pattern)
    
    # Must have at least one double or triple bond
    total_unsaturations = len(double_bonds) + len(triple_bonds)
    if total_unsaturations == 0:
        return False, "No carbon-carbon double or triple bonds found"
    
    # Count carbons to verify it's a fatty acid (typically >4 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short to be a fatty acid"
        
    # Check if the molecule has a reasonable size for a fatty acid
    # Most fatty acids have between 4-30 carbons
    if c_count > 40:
        return False, "Carbon chain too long for typical fatty acid"
        
    # Calculate the number of rotatable bonds to verify chain nature
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too rigid to be a fatty acid"
        
    # Success case - construct detailed reason
    unsaturation_details = []
    if double_bonds:
        unsaturation_details.append(f"{len(double_bonds)} double bond(s)")
    if triple_bonds:
        unsaturation_details.append(f"{len(triple_bonds)} triple bond(s)")
    
    reason = f"Contains carboxylic acid group and {' and '.join(unsaturation_details)} "
    reason += f"with {c_count} carbons in chain"
    
    return True, reason