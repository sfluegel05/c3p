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
    
    # Look for methyl ester pattern: O=C(OC)
    methyl_ester_pattern = Chem.MolFromSmarts("O=C(OC)")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"
    
    # Check for a long carbon chain (indicative of fatty acids)
    carbon_chain_pattern = Chem.MolFromSmarts("C-C-C-C-C")  # Minimal carbon chain pattern
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No carbon chain of minimum length (5) found"
    
    # Count carbon atoms - typically fatty acids have more than 8 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 9:
        return False, f"Too few carbons for a typical fatty acid methyl ester (found {c_count})"

    return True, "Contains methyl ester group and appropriate carbon chain structure for a FAME"