"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate

Definition:
  An sn-glycero-3-phosphate compound having a single unspecified acyl group
  at either position 1 or position 2. That means exactly one of the two
  hydroxyl groups (sn-1 or sn-2) in glycerol is acylated (via an ester bond)
  while the phosphate group on sn-3 is retained as an ester and no additional
  heteroatoms (e.g. nitrogen) are present.
  
Approach:
  1. Parse the SMILES.
  2. Require exactly one phosphorus atom.
  3. Reject any molecule containing nitrogen.
  4. Count acyl ester groups using the SMARTS "OC(=O)[#6]"; there must be exactly one.
  5. Look for a glycerol–phosphate backbone that has been acylated exactly at either sn‑1 or sn‑2.
     We do that with two SMARTS patterns (one for acylation at sn‑1 and one for sn‑2):
       • sn1: "[C:1](OC(=O)[C:5])[C:2]([O:4])[C:3]OP(=O)(O)O"
       • sn2: "[C:1]([O:4])[C:2](OC(=O)[C:5])[C:3]OP(=O)(O)O"
     In these patterns the free hydroxyl oxygen is explicitly tagged with mapping number 4.
  6. If a backbone match is found, ensure that the free OH (atom mapped as "4") is free,
     i.e. that in the target molecule it is bonded only to the corresponding glycerol carbon.
     
If all tests pass, the molecule is classified as monoacyl-sn-glycerol 3-phosphate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple indicating if the molecule was classified as a monoacyl-sn-glycerol 
                     3-phosphate and a reason.
    """
    # 1. Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 2. Require exactly one phosphorus atom.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, f"Molecule must contain exactly one phosphorus atom (found {len(phosphorus_atoms)})"
    
    # 3. Reject molecules containing nitrogen.
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if nitrogen_atoms:
        return False, "Molecule contains nitrogen atoms, indicating an alternative headgroup"
    
    # 4. Count acyl ester groups.
    # Define an acyl ester as an oxygen directly connected to a carbonyl carbon: OC(=O)[#6]
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    if acyl_pattern is None:
        return False, "Invalid SMARTS pattern for acyl ester"
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Molecule must contain exactly one acyl ester group (found {len(acyl_matches)})"
    
    # 5. Define SMARTS for glycerol-phosphate backbone.
    # The free hydroxyl oxygen is tagged with mapping number 4.
    sn1_smarts = "[C:1](OC(=O)[C:5])[C:2]([O:4])[C:3]OP(=O)(O)O"  # acyl at sn-1; free OH on C2
    sn2_smarts = "[C:1]([O:4])[C:2](OC(=O)[C:5])[C:3]OP(=O)(O)O"  # acyl at sn-2; free OH on C1
    patt_sn1 = Chem.MolFromSmarts(sn1_smarts)
    patt_sn2 = Chem.MolFromSmarts(sn2_smarts)
    
    if patt_sn1 is None or patt_sn2 is None:
        return False, "Invalid SMARTS pattern for glycerol-phosphate backbone"
    
    matches_sn1 = mol.GetSubstructMatches(patt_sn1)
    matches_sn2 = mol.GetSubstructMatches(patt_sn2)
    total_backbone_matches = len(matches_sn1) + len(matches_sn2)
    
    if total_backbone_matches == 0:
        return False, "No glycerol-phosphate backbone with appropriate acyl substitution found"
    if total_backbone_matches > 1:
        return False, f"Multiple ({total_backbone_matches}) backbone matches found; expected exactly one"
        
    # Identify which pattern (sn1 or sn2) provided the match.
    backbone_match = matches_sn1[0] if matches_sn1 else matches_sn2[0]
    # The backbone_match is a tuple of target atom indices corresponding to the query atoms.
    # However, to locate the atom corresponding to mapping number 4 (free hydroxyl oxygen),
    # we need to know which query index held that mapping.
    # We'll search the query (patt_sn1 if sn1 matched, else patt_sn2) for an atom with the property "molAtomMapNumber" equal to "4".
    query = patt_sn1 if matches_sn1 else patt_sn2
    free_oh_query_idx = None
    for atom in query.GetAtoms():
        if atom.HasProp("molAtomMapNumber"):
            if atom.GetProp("molAtomMapNumber") == "4":
                free_oh_query_idx = atom.GetIdx()
                break
    if free_oh_query_idx is None:
        return False, "Failed to identify free hydroxyl atom in the backbone pattern"
    
    # Now get the corresponding atom from the target molecule.
    free_oh_atom_idx = backbone_match[free_oh_query_idx]
    free_oh_atom = mol.GetAtomWithIdx(free_oh_atom_idx)
    # Check that the free OH oxygen (atomic num 8) is only bound to one heavy atom.
    heavy_neighbor_count = sum(1 for nbr in free_oh_atom.GetNeighbors() if nbr.GetAtomicNum() > 1)
    if heavy_neighbor_count != 1:
        return False, "The free hydroxyl oxygen is over-substituted; expected a free OH group"
    
    return True, "Molecule has a glycerol-phosphate backbone with exactly one acyl ester group at sn-1 or sn-2"


# For testing purposes, one may run:
if __name__ == "__main__":
    # List of example SMILES strings for monoacyl-sn-glycerol 3-phosphate molecules.
    test_smiles_list = [
        "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O",   # 1-icosanoyl-sn-glycero-3-phosphate
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCC)(O)(O)=O",      # 1-nonadecanoyl-glycero-3-phosphate
        "P(OCC(O)COC(=O)CCCCCCCCCCCCCCCCC)(O)(O)=O",      # Example: 1-stearoyl-sn-glycero-3-phosphate
    ]
    for smi in test_smiles_list:
        res, msg = is_monoacyl_sn_glycerol_3_phosphate(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {msg}\n")