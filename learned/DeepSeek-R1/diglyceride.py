"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is glycerol with two acylated hydroxy groups and one remaining group (H or alkyl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find glycerol backbone (three carbons in a row with four bonds each)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    
    # Get all possible glycerol backbones (take first match)
    glycerol_atoms = list(glycerol_matches[0])
    
    # Check for exactly two ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2]-[CX3]=[OX1]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"
    
    # Verify that both esters are attached to the glycerol backbone
    ester_oxygens = set(match[0] for match in ester_matches)
    glycerol_oxygens = []
    for atom_idx in glycerol_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() in ester_oxygens:
                glycerol_oxygens.append(neighbor.GetIdx())
    
    if len(glycerol_oxygens) != 2:
        return False, "Ester groups not attached to glycerol backbone"
    
    # Check remaining oxygen on glycerol is hydroxyl or ether
    all_glycerol_O = []
    for atom_idx in glycerol_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                all_glycerol_O.append(neighbor.GetIdx())
    
    remaining_O = [o for o in all_glycerol_O if o not in glycerol_oxygens]
    if len(remaining_O) != 1:
        return False, f"Expected 1 remaining oxygen on glycerol, found {len(remaining_O)}"
    
    remaining_o = mol.GetAtomWithIdx(remaining_O[0])
    
    # Check for hydroxyl (has H) or ether (O connected to non-carbonyl carbon)
    is_hydroxyl = remaining_o.GetTotalNumHs() > 0
    is_ether = False
    for neighbor in remaining_o.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            # Check if neighbor is not part of a carbonyl
            has_carbonyl = False
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(neighbor)
                    if other.GetAtomicNum() == 8:
                        has_carbonyl = True
                        break
            if not has_carbonyl and neighbor.GetIdx() not in glycerol_atoms:
                is_ether = True
                break
    
    if not (is_hydroxyl or is_ether):
        return False, "Remaining oxygen is neither hydroxyl nor ether"
    
    return True, "Contains glycerol backbone with two ester groups and one hydroxyl/ether group"