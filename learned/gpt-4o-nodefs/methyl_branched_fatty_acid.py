"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.

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

    # Look for carboxylic acid or carboxylate group -COOH or -COO- (smarts pattern: C(=O)[OH,[O-]])
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH,[O-]]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid or carboxylate group found"
        
    # Look for methyl branches (smarts pattern for methyl branching)
    branch_pattern = Chem.MolFromSmarts("[CH3][CH](C)")
    branch_matches = mol.GetSubstructMatches(branch_pattern)
    if not branch_matches:
        return False, "No methyl branching found"
    
    # Consider aliphatic chain along with methyl branching
    chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() == 0)
    if chain_length < 6:
        return False, "Not enough carbon atoms for a fatty acid backbone"
    
    return True, "Contains carboxylic acid or carboxylate group and methyl branching in chain"