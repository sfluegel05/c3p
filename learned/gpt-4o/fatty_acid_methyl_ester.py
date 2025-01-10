"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is characterized by an ester linkage from methanol 
    and an aliphatic carbon chain similar to fatty acids, potentially with unsaturations or branching.

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

    # Check for methyl ester group (methanol-derived ester link)
    ester_methyl_pattern = Chem.MolFromSmarts("O=C(OC)")
    if not mol.HasSubstructMatch(ester_methyl_pattern):
        return False, "Methyl ester group missing"

    # Look for long carbon chains, accounting for branching or double bonds
    long_chain_pattern = Chem.MolFromSmarts("C(C)(C)(C)C")  # Detect at least a 4-carbon chain allowing branches
    unsat_chain_pattern = Chem.MolFromSmarts("C=CC")        # Detect carbon double bonds in chains
    
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    unsat_chain_matches = mol.GetSubstructMatches(unsat_chain_pattern)
    
    # A valid fatty acid methyl ester should have sufficient carbon length
    # and optionally unsaturation or branches
    if len(long_chain_matches) < 3 and len(unsat_chain_matches) < 1:  # Adjusted to consider >3 carbons more flexibly
        return False, "Missing sufficient long carbon chain or unsaturation for fatty acid"

    return True, "Structure matches fatty acid methyl ester requirements"