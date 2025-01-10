"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a fatty acid containing only methyl branches.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should have just one
    carboxylic_matches = len(mol.GetSubstructMatches(carboxylic_pattern))
    if carboxylic_matches != 1:
        return False, f"Invalid number of carboxylic acid groups ({carboxylic_matches})"
    
    # Check for aromatic character
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Find the carbon attached to COOH group
    acid_carbon_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    matches = mol.GetSubstructMatches(acid_carbon_pattern)
    if not matches:
        return False, "Cannot identify acid carbon"
    
    # Analyze branching patterns
    branch_patterns = {
        "ethyl": "[CH2][CH2][CH2,CH3]",  # Ethyl or longer branches
        "propyl": "[CH2][CH2][CH2][CH2,CH3]",  # Propyl or longer branches
        "methyl": "[CH3][CX4]",  # Methyl branches
        "isopropyl": "[CH3][CH]([CH3])",  # Isopropyl groups
    }
    
    # Count different types of branches
    branch_counts = {}
    for name, pattern in branch_patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        matches = len(mol.GetSubstructMatches(pat))
        branch_counts[name] = matches

    # Check for non-methyl branches
    if branch_counts["ethyl"] > 1 or branch_counts["propyl"] > 0 or branch_counts["isopropyl"] > 0:
        return False, "Contains non-methyl branches"
    
    # Must have at least one methyl branch
    if branch_counts["methyl"] == 0:
        return False, "No methyl branches found"
    
    # Count carbons in longest chain
    longest_chain = rdMolDescriptors.CalcMolFormula(mol).count('C')
    if longest_chain < 4:
        return False, "Carbon chain too short"
    if longest_chain > 30:
        return False, "Carbon chain too long for typical fatty acid"

    # Check heteroatom content
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    
    # Allow only C, H, O (and limited N)
    if any(symbol not in ['C', 'H', 'O', 'N'] for symbol in atom_counts.keys()):
        return False, "Contains unexpected heteroatoms"
    
    # Check nitrogen content
    if atom_counts.get('N', 0) > 1:
        return False, "Too many nitrogen atoms"

    # Check for specific structural features that would disqualify
    disqualifying_patterns = [
        "[N+]",  # Quaternary nitrogen
        "[S,P,B,Si]",  # Other heteroatoms
        "C(=O)OC(=O)",  # Anhydrides
        "[OH]C(=O)[OH]"  # Carbonic acid
    ]
    
    for pattern in disqualifying_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Contains disqualifying structural feature: {pattern}"

    return True, "Contains carboxylic acid group with methyl branching"