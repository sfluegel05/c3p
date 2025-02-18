"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine (CHEBI:89834)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    Must have a glycerol backbone with two acyl groups at positions 1 and 2, and a phosphoserine group at position 3.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule matches criteria, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # 1. Find glycerol backbone with sn-3 configuration (C2 is the central carbon)
    # Pattern: [CH2]-[C@H](-O-...)-[CH2]-O-...
    glycerol_pattern = Chem.MolFromSmarts("[CH2][C@H]([OX2])[CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No sn-glycerol backbone"
    
    # Get the central carbon (C2) of the glycerol
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "Glycerol pattern not found"
    c2_idx = matches[0][1]  # Central carbon is the second atom in the match
    
    # 2. Check substituents on C2: should have two ester oxygens and one phosphate oxygen
    c2 = mol.GetAtomWithIdx(c2_idx)
    ester_count = 0
    phosphate_oxygen = None
    
    for neighbor in c2.GetNeighbors():
        if neighbor.GetSymbol() == "O":
            # Check if the oxygen is part of an ester group (O-C=O)
            for bond in neighbor.GetBonds():
                if bond.GetBeginAtomIdx() == neighbor.GetIdx() or bond.GetEndAtomIdx() == neighbor.GetIdx():
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetSymbol() == "C" and any(a.GetSymbol() == "O" for a in other_atom.GetNeighbors() if a.GetBondType() == Chem.BondType.DOUBLE):
                        ester_count += 1
            # Check if the oxygen is connected to phosphorus
            for bond in neighbor.GetBonds():
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetSymbol() == "P":
                    phosphate_oxygen = neighbor
                    break
    
    if ester_count != 2:
        return False, f"Expected 2 ester groups on C2, found {ester_count}"
    if not phosphate_oxygen:
        return False, "No phosphate group attached to C2"
    
    # 3. Verify phosphoserine group: follow the phosphate oxygen to serine
    # Phosphoserine should have NH2 and COOH groups connected via the phosphate's oxygen chain
    p_atom = phosphate_oxygen.GetNeighbors()[0]
    # Traverse from P to find serine moiety: P-O-C-C(N)-C(=O)O
    # Pattern for serine part: C connected to NH2 and COOH
    serine_pattern = Chem.MolFromSmarts("[NH2]-C(-C(=O)[OH])-C-O-P")
    if not mol.HasSubstructMatch(serine_pattern):
        # Alternative pattern considering possible charges
        serine_pattern_alt = Chem.MolFromSmarts("[NH2]-C(-C(=O)[O-])-C-O-P")
        if not mol.HasSubstructMatch(serine_pattern_alt):
            return False, "Phosphoserine group not found"
    
    # 4. Check acyl chains at positions 1 and 2 (minimum length)
    # Get ester groups attached to C2's oxygens
    # Assuming the two ester groups are the acyl chains
    # Check for at least 4 carbons in each chain (arbitrary threshold for "long" chain)
    ester_groups = []
    for neighbor in c2.GetNeighbors():
        if neighbor.GetSymbol() == "O":
            for bond in neighbor.GetBonds():
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetSymbol() == "C" and any(a.GetSymbol() == "O" for a in other_atom.GetNeighbors() if a.GetBondType() == Chem.BondType.DOUBLE):
                    ester_groups.append(other_atom)
    
    if len(ester_groups) < 2:
        return False, "Not enough ester groups"
    
    # Check chain length for each ester
    min_chain_length = 4  # Adjust based on requirements
    for ester in ester_groups[:2]:  # First two esters (positions 1 and 2)
        chain = []
        current = ester
        prev = ester.GetNeighbors()[0]  # O atom
        while True:
            next_atoms = [a for a in current.GetNeighbors() if a != prev and a.GetSymbol() in ["C", "H"]]
            if not next_atoms or current.GetSymbol() != "C":
                break
            chain.append(current)
            prev = current
            current = next_atoms[0]
        if len(chain) < min_chain_length:
            return False, f"Acyl chain too short: {len(chain)} carbons"
    
    return True, "sn-glycerol with two acyl chains and phosphoserine at position 3"