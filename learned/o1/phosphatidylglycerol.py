"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol consists of a glycerol backbone with two fatty acid chains attached via ester bonds,
    a phosphate group attached to the third carbon, and another glycerol attached to the phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern with ester bonds on C1 and C2, phosphate on C3
    glycerol_pattern = Chem.MolFromSmarts("""
    [$([C@@H]1([O][P](=O)([O])[O][C@@H]([C@@H]([O])CO)O)[O][C](=O)C)] 
    -[$([O][C](=O)C)] 
    -[C@@H]([O])CO
    """)
    if not glycerol_pattern:
        return False, "Error in glycerol SMARTS pattern"

    # Check for the glycerol backbone with attached groups
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No phosphatidylglycerol core structure found"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("OC(=O)[C;H2,C;H1,C;H0]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Look for phosphate group attached to glycerol backbone
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group attached to glycerol backbone"

    # Look for second glycerol attached to phosphate group
    glycerol_phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(OC[C@@H](CO)O)")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol attached to phosphate group"

    # Optional: Check for long carbon chains in fatty acids
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)OCCCCC")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No long-chain fatty acids found"

    # All structural features matched
    return True, "Molecule is a phosphatidylglycerol with correct structural features"