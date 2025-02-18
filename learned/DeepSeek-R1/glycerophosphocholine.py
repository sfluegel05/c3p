"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: CHEBI:17585 glycerophosphocholine
"""
from rdkit import Chem
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

    # Check glycerol backbone connected to phosphate
    # Glycerol pattern: C-O-C-O-C where middle C is connected to phosphate
    # SMARTS for glycerol part connected to phosphate: [CH2]-O-[CH](OP(=O)([O-])OCC[N+])-O-[CH2]
    # More precise pattern to ensure correct connectivity
    glycerol_phosphate = Chem.MolFromSmarts("[CH2]-[O]-[CH](-[O]-P(=O)([O-])-[O]-C-C-[N+](C)(C)C)-[O]-[CH2]")
    if not mol.HasSubstructMatch(glycerol_phosphate):
        return False, "Glycerol-phosphate backbone not found"

    # Check for two ester or ether-linked fatty acid chains
    ester_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-C(=O)")  # Ester group
    ether_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")  # Ether group (no carbonyl)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    total_fatty = len(ester_matches) + len(ether_matches)
    
    # Subtract the phosphate oxygen (part of phosphocholine) from ether matches
    # Phosphate has three oxygens: two in P=O and one connecting to glycerol
    # The ether oxygen in the phosphate group is part of the backbone, not a fatty chain
    # Adjust count by subtracting the phosphate's ether oxygen
    adjusted_total = total_fatty - 1  # subtract the phosphate oxygen
    if adjusted_total < 2:
        return False, f"Found {adjusted_total} fatty acid chains (need 2)"

    # Check molecular weight (adjust threshold as needed)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:  # Lowered threshold to accommodate shorter chains
        return False, f"Molecular weight {mol_wt:.1f} too low"

    return True, "Glycerol backbone with phosphocholine and two fatty acid chains"