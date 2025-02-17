"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids

Definition:
    Lysophosphatidic acids (LPA) are monoacylglycerol phosphates derived from phosphatidic
    acid by removal of one of the acyl groups. They should contain a phosphate group directly
    linked to a glycerol backbone and exactly one acyl ester bond (that attaches a fatty acid).
    
Improved heuristic criteria:
    1. The molecule must contain a phosphate group with the pattern P(=O)(O)(O).
    2. At least one of the oxygen atoms of the phosphate must be directly bonded to a carbon
       that belongs to a glycerol backbone.
    3. The molecule must contain a glycerol fragment. We use two alternative SMARTS patterns to be flexible.
    4. The molecule must contain exactly ONE ester bond (substructure “OC(=O)”) for attachment of one acyl chain.
       We require that the “ester oxygen” is part of the glycerol backbone.
    5. Its molecular weight must be above ~250 Da (to filter out small molecules).
Note: This is only a heuristic and will not be 100% specific.
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
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 1. Check for the phosphate group string.
    # Use a SMARTS pattern that matches a phosphate: P(=O)(O)(O)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    phos_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phos_matches:
        return False, "Phosphate group not found."
    
    # 2. Check that at least one oxygen bound to the phosphate connects to a glycerol backbone.
    # We look for the glycerol fragment. Because stereochemistry may be present or not,
    # we use two patterns.
    # First pattern: a simplified glycerol backbone: O C C(O) C O  (ignoring stereochemistry)
    glycerol_pattern1 = Chem.MolFromSmarts("OCC(O)CO")
    # Alternative pattern (common when chiral tags are given)
    glycerol_pattern2 = Chem.MolFromSmarts("OC[C@H](O)CO")
    
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern1)
    if not glycerol_matches:
        glycerol_matches = mol.GetSubstructMatches(glycerol_pattern2)
    if not glycerol_matches:
        return False, "Glycerol backbone not found."
    # For later use, collect all atom indices in any glycerol match (union over each match)
    glycerol_atoms = set()
    for match in glycerol_matches:
        glycerol_atoms.update(match)
    
    # 3. Verify that the phosphate group is attached to the glycerol fragment.
    # For each phosphorus-containing match, check if any oxygen neighbor is bonded to a carbon that’s in glycerol.
    valid_phosphate_found = False
    for match in phos_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 15:  # phosphorus
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                        # if one of oxygen's neighbors (except the phosphorus) is in the glycerol fragment, we consider it attached.
                        for nbr2 in nbr.GetNeighbors():
                            if nbr2.GetIdx() in glycerol_atoms and nbr2.GetAtomicNum() == 6:
                                valid_phosphate_found = True
                                break
                    if valid_phosphate_found:
                        break
            if valid_phosphate_found:
                break
        if valid_phosphate_found:
            break
    if not valid_phosphate_found:
        return False, "Phosphate group not directly linked to a glycerol backbone."

    # 4. Identify the ester bond: the fatty acid attachment.
    # We look for the ester substructure "OC(=O)".
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    # We now count only those ester matches where the oxygen (first atom in the match) is part of the glycerol backbone.
    ester_count = 0
    for match in ester_matches:
        oxygen_idx = match[0]
        if oxygen_idx in glycerol_atoms:
            ester_count += 1
    if ester_count != 1:
        return False, f"Found {ester_count} ester bond(s) attached to glycerol; expected exactly 1 for monoacyl LPA."

    # 5. Check that the molecular weight is within an expected range for an LPA molecule.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low to be a lysophosphatidic acid."
    
    return True, "Molecule contains a phosphate group directly linked to a glycerol backbone with exactly one acyl ester bond."

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples.
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(O)(O)=O"  # PA(17:1(9Z)/0:0)
    result, reason = is_lysophosphatidic_acids(test_smiles)
    print(result, reason)