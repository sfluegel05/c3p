"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids

Definition:
    Lysophosphatidic acids (LPA) are monoacylglycerol phosphates obtained by removal of one
    acyl group from phosphatidic acid. In our heuristic, an LPA should:
      1. Contain a phosphate group with the pattern P(=O)(O)(O).
      2. Have the phosphate directly attached to a glycerol backbone.
      3. Contain exactly one fatty acid ester bond (the “O–C(=O)” bond) attached to the glycerol.
      4. Have a molecular weight above ~250 Da.
      5. Not contain additional nitrogen atoms (which often indicate phosphocholines or similar lipids).
Note: This is a heuristic method and may not be 100% specific.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid (LPA) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as an LPA, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Reject molecules that contain nitrogen atoms,
    # since LPAs should not have nitrogen (they lack choline/ethanolamine groups).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            return False, "Molecule contains nitrogen atoms, likely a phosphocholine/phosphoethanolamine."
    
    # 1. Check for the phosphate group.
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    phos_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phos_matches:
        return False, "Phosphate group not found."
    
    # 2. Find glycerol backbone fragments.
    # We use two patterns to allow for different representations (with/without explicit stereochemistry).
    glycerol_pattern1 = Chem.MolFromSmarts("OCC(O)CO")
    glycerol_pattern2 = Chem.MolFromSmarts("OC[C@H](O)CO")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern1)
    if not glycerol_matches:
        glycerol_matches = mol.GetSubstructMatches(glycerol_pattern2)
    if not glycerol_matches:
        return False, "Glycerol backbone not found."
    
    # Combine all atom indices that are part of any glycerol match.
    glycerol_atoms = set()
    for match in glycerol_matches:
        glycerol_atoms.update(match)
    
    # 3. Verify that the phosphate group is attached to the glycerol backbone.
    # For each phosphate match, look at the phosphorus atom’s oxygen neighbors.
    found_linkage = False
    for match in phos_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 15:  # phosphorus atom
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                        # Check if oxygen's neighbors (except the P) include a carbon in the glycerol fragment.
                        for nbr2 in nbr.GetNeighbors():
                            if nbr2.GetIdx() != atom_idx and nbr2.GetAtomicNum() == 6:
                                if nbr2.GetIdx() in glycerol_atoms:
                                    found_linkage = True
                                    break
                        if found_linkage:
                            break
                if found_linkage:
                    break
        if found_linkage:
            break
    if not found_linkage:
        return False, "Phosphate group not directly linked to glycerol backbone."
    
    # 4. Identify the ester bond (attachment of one fatty acyl chain).
    # We search for the ester substructure "OC(=O)".
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ester_count = 0
    for match in ester_matches:
        # The first atom (an oxygen) in the match should belong to the glycerol backbone.
        oxygen_idx = match[0]
        if oxygen_idx in glycerol_atoms:
            ester_count += 1
    if ester_count != 1:
        return False, f"Found {ester_count} ester bond(s) attached to glycerol; expected exactly 1 for monoacyl LPA."
    
    # 5. Check the molecular weight (filter out very small molecules).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low to be a lysophosphatidic acid."
    
    return True, "Molecule contains a phosphate group directly linked to a glycerol backbone with exactly one acyl ester bond."

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES string.
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(O)(O)=O"  # PA(17:1(9Z)/0:0)
    result, reason = is_lysophosphatidic_acids(test_smiles)
    print(result, reason)