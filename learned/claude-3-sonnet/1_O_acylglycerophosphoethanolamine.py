"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    These compounds have:
    - A glycerol backbone
    - An acyl group (ester) at sn-1 position
    - A hydroxyl at sn-2 position
    - A phosphoethanolamine group at sn-3 position
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure pattern that matches both neutral and zwitterionic forms
    # Includes glycerol backbone with ester at sn-1, OH at sn-2, and phosphoethanolamine at sn-3
    core_pattern = Chem.MolFromSmarts("""
        [C,c:1](=[O:2])[O:3][CH2:4][CH:5]([OH1,O:6])[CH2:7][O:8]P(=[O:9])([O,OH:10])[O:11][CH2:12][CH2:13][NH2,NH3+:14]
    """)
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing core glycerophosphoethanolamine structure"

    # Count key features to ensure we have the right structure
    # Count phosphorus atoms - should be exactly 1
    p_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#15]")))
    if p_count != 1:
        return False, f"Should have exactly 1 phosphorus, found {p_count}"

    # Count ester groups - should be exactly 1
    ester_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[C,c](=O)[O][CH2,CH]")))
    if ester_count != 1:
        return False, f"Should have exactly 1 ester group, found {ester_count}"

    # Count nitrogen atoms - should be exactly 1
    n_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#7]")))
    if n_count != 1:
        return False, f"Should have exactly 1 nitrogen, found {n_count}"

    # Verify it's not a phosphocholine (no N(CH3)3+)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NX4+](C)(C)(C)")):
        return False, "Contains quaternary amine (phosphocholine)"

    # Check for minimum acyl chain length (at least 4 carbons)
    acyl_chain = Chem.MolFromSmarts("[C,c](=O)[O][CH2][CH]([OH1,O])[CH2][O]P")
    matches = mol.GetSubstructMatches(acyl_chain)
    if not matches:
        return False, "Could not verify acyl chain position"
    
    # Get the first carbon of the acyl group and count connected carbons
    for match in matches:
        acyl_c = mol.GetAtomWithIdx(match[0])
        chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if chain_length < 7:  # Minimum total carbons needed (3 glycerol + 4 acyl chain)
            return False, "Acyl chain too short"

    return True, "Contains glycerol backbone with phosphoethanolamine group and single acyl substituent at sn-1"