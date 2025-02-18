"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: CHEBI:175825 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    Must have a glycerol backbone with:
    - 1-O-acyl group (ester at position 1)
    - Phosphoethanolamine group at position 3
    - Only one acyl chain (lyso structure)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphoethanolamine group (O-P-OCCN)
    pe_pattern = Chem.MolFromSmarts("[O][P](=O)([O][C][C][N])[O][C]")
    if not mol.HasSubstructMatch(pe_pattern):
        return False, "Missing phosphoethanolamine group"

    # Find glycerol backbone connected to phosphate
    # Glycerol pattern: C-O-C-O-P with ester at first carbon
    core_pattern = Chem.MolFromSmarts("[CH2]-[CH](-O)-[CH2]-O-P(=O)(OCCN)")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing glycerol-phosphate core"

    # Verify ester group at position 1 (O-C=O attached to first glycerol carbon)
    ester_pattern = Chem.MolFromSmarts("[CH2][OX2][C](=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check acyl chain length (minimum 10 carbons including carbonyl)
    carbonyl = Chem.MolFromSmarts("[CX3]=[OX1]")
    for match in mol.GetSubstructMatches(carbonyl):
        atom = mol.GetAtomWithIdx(match[0])
        chain_length = 0
        # Traverse connected carbons excluding the glycerol oxygen
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != match[1]:
                chain_length = 1 + rdMolDescriptors.CalcNumAdjacentHeavyAtoms(neighbor)
                break
        if chain_length < 10:
            return False, f"Acyl chain too short ({chain_length} carbons)"

    # Check no other ester groups on glycerol
    all_esters = mol.GetSubstructMatches(Chem.MolFromSmarts("[O][C](=O)"))
    if len(all_esters) > 1:
        return False, "Multiple acyl groups detected"

    return True, "1-O-acylglycerophosphoethanolamine structure confirmed"