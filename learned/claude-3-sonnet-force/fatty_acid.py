"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:36976 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is any aliphatic monocarboxylic acid derived from or contained in
    esterified form in an animal or vegetable fat, oil or wax. Natural fatty acids
    commonly have a chain of 4 to 28 carbons (usually unbranched and even-numbered),
    which may be saturated or unsaturated. By extension, the term is sometimes used
    to embrace all acyclic aliphatic carboxylic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for aliphatic segments (potentially interrupted by double bonds or functional groups)
    aliphatic_pattern = Chem.MolFromSmarts("[C;H3,H2,H1]~[C;H3,H2,H1]")
    aliphatic_segments = mol.GetSubstructMatches(aliphatic_pattern)
    
    total_aliphatic_carbons = 0
    for segment in aliphatic_segments:
        # Expand segment in both directions to find full aliphatic chain
        start, end = segment[0], segment[-1]
        chain = AllChem.FindAllPathsOfLengthN(mol, 1000, start, end)[0]
        total_aliphatic_carbons += len(chain)
    
    if total_aliphatic_carbons < 4 or total_aliphatic_carbons > 28:
        return False, f"Total aliphatic carbon count of {total_aliphatic_carbons} is outside typical fatty acid range (4-28)"
    
    # Check for common functional groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    epoxide_pattern = Chem.MolFromSmarts("[O;r5]")
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    has_halogen = mol.HasSubstructMatch(halogen_pattern)
    
    # Check molecular weight and size
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 800:
        return False, "Molecular weight outside typical range for fatty acids"
    
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 10 or n_atoms > 60:
        return False, "Number of atoms outside typical range for fatty acids"
    
    functional_groups = []
    if has_hydroxy:
        functional_groups.append("hydroxy")
    if has_epoxide:
        functional_groups.append("epoxide")
    if has_halogen:
        functional_groups.append("halogen")
    
    functional_group_str = ", ".join(functional_groups) if functional_groups else "none"
    
    return True, f"Contains carboxylic acid group and aliphatic segment(s) with {total_aliphatic_carbons} carbons. Functional groups: {functional_group_str}"