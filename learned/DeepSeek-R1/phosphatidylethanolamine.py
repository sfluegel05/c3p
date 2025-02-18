"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine has a glycerol backbone with two fatty acid esters and a phosphoethanolamine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphoethanolamine group: P connected to O-C-C-N (with possible N-methylation)
    pe_pattern = Chem.MolFromSmarts("[PX4](=O)(-[OX2][CX4][CX4][NX3;H2,H1,H0])")
    if not mol.HasSubstructMatch(pe_pattern):
        return False, "No phosphoethanolamine group detected"

    # Check for two ester groups (O-C=O)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    esters = mol.GetSubstructMatches(ester_pattern)
    if len(esters) < 2:
        return False, f"Found {len(esters)} ester groups, need at least 2"

    # Verify glycerol backbone connectivity
    # Find phosphorus atoms in phosphoethanolamine groups
    pe_matches = mol.GetSubstructMatches(pe_pattern)
    for match in pe_matches:
        p_atom = mol.GetAtomWithIdx(match[0])
        
        # Find oxygen connecting to glycerol
        glycerol_o = None
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O' and not neighbor.GetTotalNumHs():
                glycerol_o = neighbor
                break
        if not glycerol_o:
            continue
        
        # Find connected carbon (glycerol backbone)
        glycerol_c = None
        for neighbor in glycerol_o.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                glycerol_c = neighbor
                break
        if not glycerol_c:
            continue
        
        # Check for two ester-connected oxygens on adjacent carbons
        ester_count = 0
        for neighbor in glycerol_c.GetNeighbors():
            if neighbor.GetSymbol() == 'O':
                # Check if oxygen is part of an ester group
                for o_neighbor in neighbor.GetNeighbors():
                    if o_neighbor.GetSymbol() == 'C' and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in o_neighbor.GetBonds()):
                        ester_count += 1
                        break
        
        # Check adjacent carbons in glycerol chain
        adjacent_ester_count = 0
        for neighbor in glycerol_c.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                for c_neighbor in neighbor.GetNeighbors():
                    if c_neighbor.GetSymbol() == 'O':
                        for o_neighbor in c_neighbor.GetNeighbors():
                            if o_neighbor.GetSymbol() == 'C' and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in o_neighbor.GetBonds()):
                                adjacent_ester_count += 1
        
        if ester_count + adjacent_ester_count >= 2:
            return True, "Contains glycerol backbone with two ester groups and phosphoethanolamine"

    return False, "Phosphoethanolamine group not properly connected to glycerol with two esters"