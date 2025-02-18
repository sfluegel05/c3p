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

    # Improved phosphoethanolamine pattern: P connected to O-C-C-N (with possible substituents on N)
    pe_pattern = Chem.MolFromSmarts("[PX4](=O)(O-[CX4][CX4][NX3;H2,H1,H0])")
    pe_matches = mol.GetSubstructMatches(pe_pattern)
    if not pe_matches:
        return False, "No phosphoethanolamine group detected"

    # Check each phosphoethanolamine group's connection to glycerol
    for pe_match in pe_matches:
        p_idx = pe_match[0]
        p_atom = mol.GetAtomWithIdx(p_idx)
        
        # Find oxygen connecting to glycerol (should be single bond to C)
        glycerol_o = None
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 0:
                for bond in p_atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.SINGLE and bond.GetOtherAtom(p_atom) == neighbor:
                        glycerol_o = neighbor
                        break
                if glycerol_o:
                    break
        if not glycerol_o:
            continue
        
        # Get connected carbon (C3 of glycerol)
        c3 = glycerol_o.GetNeighbors()[0]
        if c3.GetSymbol() != 'C':
            continue
        
        # Trace glycerol backbone C1-C2-C3
        # Find C2 connected to C3 (excluding the phosphate oxygen)
        c2 = None
        for neighbor in c3.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != glycerol_o.GetIdx():
                c2 = neighbor
                break
        if not c2:
            continue
        
        # Find C1 connected to C2 (excluding C3)
        c1 = None
        for neighbor in c2.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != c3.GetIdx():
                c1 = neighbor
                break
        if not c1:
            continue
        
        # Check ester groups on C1 and C2
        ester_count = 0
        
        # Check C1 for ester
        for bond in c1.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtom(c1)
                if neighbor.GetSymbol() == 'O':
                    for o_bond in neighbor.GetBonds():
                        o_neighbor = o_bond.GetOtherAtom(neighbor)
                        if o_neighbor.GetSymbol() == 'C' and any(b.GetBondType() == Chem.BondType.DOUBLE for b in o_neighbor.GetBonds()):
                            ester_count +=1
                            break
        
        # Check C2 for ester
        for bond in c2.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtom(c2)
                if neighbor.GetSymbol() == 'O':
                    for o_bond in neighbor.GetBonds():
                        o_neighbor = o_bond.GetOtherAtom(neighbor)
                        if o_neighbor.GetSymbol() == 'C' and any(b.GetBondType() == Chem.BondType.DOUBLE for b in o_neighbor.GetBonds()):
                            ester_count +=1
                            break
        
        if ester_count >= 2:
            return True, "Glycerol backbone with two esters and phosphoethanolamine group"
    
    return False, "Missing required glycerol-ester-phosphoethanolamine structure"