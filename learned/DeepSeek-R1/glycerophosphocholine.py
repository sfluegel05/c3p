"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: CHEBI:17585 glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine consists of a glycerol backbone with a phosphocholine group and two fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphocholine group (O-P-O-C-C-N+(C)(C)C)
    phosphocholine_pattern = Chem.MolFromSmarts("[O]-P(=O)([O-])-O-C-C-[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine group not found"

    # Identify the glycerol backbone (C-O-C-O-C with phosphate attachment)
    # The glycerol is typically C1-O-C2-O-C3, where C2 is connected to the phosphate
    # Look for a central carbon connected to two ether/ester oxygens and one phosphate oxygen
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[O]-[C;!$(C=O)]")  # Adjust as needed
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not identified"

    # Check for two ester/ether-linked fatty acid chains
    ester_ether = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4](=O) | [CX4]-[OX2]-[CX4]")  # Ester or ether
    matches = mol.GetSubstructMatches(ester_ether)
    if len(matches) < 2:
        return False, f"Found {len(matches)} ester/ether groups, need at least 2"

    # Verify molecular weight is consistent with fatty acids (optional)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 600:  # Adjust threshold based on typical glycerophosphocholines
        return False, f"Molecular weight {mol_wt:.1f} too low for glycerophosphocholine"

    return True, "Contains glycerol backbone with phosphocholine and two fatty acid chains"