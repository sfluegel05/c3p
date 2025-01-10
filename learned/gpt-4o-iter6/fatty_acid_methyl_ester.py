"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester (FAME) based on its SMILES string.
    A FAME is defined as a carboxylic ester obtained by the formal condensation of a fatty acid with methanol.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ester group with a methyl: O=C(OC)
    # Considering simple esters bonded to a single carboxyl,
    # and avoiding cases where methyl is deeply embedded in larger structures.
    methyl_ester_pattern = Chem.MolFromSmarts("O=C(OC);-[C!$(*=*)C!$(*=*)]-;!")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"
    
    # Check for a carbon chain indicative of fatty acids
    # Allow for unsaturation and slight branching.
    carbon_chain_flexible_pattern = Chem.MolFromSmarts("C!@[C,!$(CC)]*C-;!@[C!$(CC)][C;!$(*=*)]C")
    
    # We also check for consecutive flexible atoms sequence to avoid
    # misclassification by trivial short chains or simple carbon groups.
    chain_matches = mol.GetSubstructMatches(carbon_chain_flexible_pattern)
    if not chain_matches:
        return False, "No indicative carbon chain found"

    # Calculate carbon atoms outside of ester moiety
    # A single carbon should be methyl ester moiety, separated from the chain.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:  # Account for at least 5 carbons beyond esterification
        return False, f"Too few carbons for a fatty acid methyl ester (found {c_count})"

    return True, "Contains methyl ester group and a suitable carbon chain for a FAME"