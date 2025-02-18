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
    
    # 1. Find glycerol backbone with sn-3 configuration (central carbon with specified stereochemistry)
    # Pattern matches: [CH2]-[C@H](-O-*)-[CH2] where * represents substituents
    glycerol_pattern = Chem.MolFromSmarts("[CH2][C@H]([OX2])[CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No sn-glycerol backbone found"
    
    # Get the central carbon (C2) of the glycerol
    matches = mol.GetSubstructMatches(glycerol_pattern)
    c2_idx = matches[0][1]  # Central carbon is the second atom in the match
    
    # 2. Check substituents on C2: should have two ester oxygens and one phosphate oxygen
    c2 = mol.GetAtomWithIdx(c2_idx)
    ester_count = 0
    phosphate_found = False
    
    for neighbor in c2.GetNeighbors():
        if neighbor.GetSymbol() == "O":
            # Check if oxygen is part of an ester group (O-C(=O))
            ester_flag = False
            for bond in neighbor.GetBonds():
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetSymbol() == "C":
                    # Check if this carbon has a double bond to oxygen
                    for a in other_atom.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(other_atom.GetIdx(), a.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE and a.GetSymbol() == "O":
                            ester_count += 1
                            ester_flag = True
                            break
                    if ester_flag:
                        break
            
            # Check if oxygen is part of phosphate group
            for bond in neighbor.GetBonds():
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetSymbol() == "P":
                    phosphate_found = True
    
    if ester_count != 2:
        return False, f"Expected 2 ester groups, found {ester_count}"
    if not phosphate_found:
        return False, "No phosphate group attached to glycerol"
    
    # 3. Verify phosphoserine group using SMARTS pattern
    # Pattern matches: P-O-C-C(N)-C(=O)O with possible charges
    serine_pattern = Chem.MolFromSmarts("[P](=O)(O)OCC([NH2])C(=O)O")
    if not mol.HasSubstructMatch(serine_pattern):
        # Try alternative pattern with carboxylate form
        serine_pattern_alt = Chem.MolFromSmarts("[P](=O)(O)OCC([NH2])C(=O)[O-]")
        if not mol.HasSubstructMatch(serine_pattern_alt):
            return False, "Phosphoserine group not found"
    
    # 4. Verify acyl chain lengths (at least 4 carbons each)
    ester_groups = []
    for neighbor in c2.GetNeighbors():
        if neighbor.GetSymbol() == "O":
            for bond in neighbor.GetBonds():
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetSymbol() == "C":
                    # Check for carbonyl group
                    for a in other_atom.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(other_atom.GetIdx(), a.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE and a.GetSymbol() == "O":
                            ester_groups.append(other_atom)
                            break
    
    if len(ester_groups) < 2:
        return False, "Insufficient ester groups for acyl chains"
    
    # Count chain length using molecular weight as proxy
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:  # Typical phosphatidylserines are >600 Da
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"
    
    return True, "sn-glycerol with two acyl chains and phosphoserine at position 3"