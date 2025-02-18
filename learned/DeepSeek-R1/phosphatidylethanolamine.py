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

    # Match phosphoethanolamine group (P=O connected to O-C-C-N with possible N-methylation)
    pe_pattern = Chem.MolFromSmarts("[PX4](=O)(O-[CX4]-[CX4]-[NX3;H2,H1,H0])")
    pe_matches = mol.GetSubstructMatches(pe_pattern)
    if not pe_matches:
        return False, "No phosphoethanolamine group detected"

    # Find all ester groups in the molecule
    ester_pattern = Chem.MolFromSmarts("[OX2]-C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Check each phosphoethanolamine group's connection to glycerol backbone
    for pe_match in pe_matches:
        p_idx = pe_match[0]
        p_atom = mol.GetAtomWithIdx(p_idx)
        
        # Find oxygen connecting phosphate to glycerol
        glycerol_o = None
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 0:
                # Verify single bond between P and O
                bond = mol.GetBondBetweenAtoms(p_idx, neighbor.GetIdx())
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    glycerol_o = neighbor
                    break
        if not glycerol_o:
            continue
        
        # Get connected carbon (C3 of glycerol)
        c3 = glycerol_o.GetNeighbors()[0]
        if c3.GetSymbol() != 'C':
            continue
        
        # Find adjacent carbons (C2 and C1) to form glycerol backbone C1-C2-C3
        c2 = None
        for neighbor in c3.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != glycerol_o.GetIdx():
                c2 = neighbor
                break
        if not c2:
            continue
        
        c1 = None
        for neighbor in c2.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != c3.GetIdx():
                c1 = neighbor
                break
        if not c1:
            continue
        
        # Collect all glycerol carbons
        glycerol_carbons = {c1.GetIdx(), c2.GetIdx(), c3.GetIdx()}
        
        # Count esters attached to glycerol carbons
        glycerol_ester_count = 0
        for ester_match in ester_matches:
            ester_o_idx = ester_match[0]
            ester_o = mol.GetAtomWithIdx(ester_o_idx)
            # Check if ester oxygen is connected to a glycerol carbon
            for neighbor in ester_o.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() in glycerol_carbons:
                    glycerol_ester_count += 1
                    break
        
        if glycerol_ester_count >= 2:
            return True, "Glycerol backbone with two esters and phosphoethanolamine group"
    
    return False, "Missing required glycerol-ester-phosphoethanolamine structure"