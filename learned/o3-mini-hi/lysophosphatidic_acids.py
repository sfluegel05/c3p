"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids

Definition:
    Lysophosphatidic acids (LPA) are monoacylglycerol phosphates obtained by removal of one
    acyl group from phosphatidic acid. In our heuristic an LPA should:
      1. Contain a phosphate group with the pattern P(=O)(O)(O).
      2. Have the phosphate directly attached to a glycerol backbone.
      3. Contain exactly one fatty acid ester bond (i.e. an “O–C(=O)” linkage) attached to the glycerol.
      4. Have a molecular weight above ~250 Da.
      5. Not contain nitrogen atoms.
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
        str: Explanation text for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 0. Reject molecules that contain nitrogen atoms. (LPAs normally lack nitrogen.)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            return False, "Molecule contains nitrogen atoms; likely not an LPA."
    
    # 1. Identify the phosphate group.
    # The phosphate pattern here is chosen as P(=O)(O)(O) (it can appear deprotonated as well).
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    phos_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phos_matches:
        return False, "Phosphate group not found."

    # 2. Identify glycerol backbone fragments.
    # We try two slightly different SMARTS patterns to capture different representations.
    glycerol_pattern1 = Chem.MolFromSmarts("OCC(O)CO")
    glycerol_pattern2 = Chem.MolFromSmarts("OC[C@H](O)CO")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern1)
    if not glycerol_matches:
        glycerol_matches = mol.GetSubstructMatches(glycerol_pattern2)
    if not glycerol_matches:
        return False, "Glycerol backbone not found."
    
    # 3. Identify ester bonds.
    # We use the SMARTS for an ester linkage "OC(=O)".
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    # To ease later lookup, we will build a list of ester oxygen indices with their full match.
    ester_oxygens = [match[0] for match in ester_matches]  # first atom of match is the linking oxygen

    # 4. Now analyze each glycerol backbone match to see if it qualifies.
    # We expect one of the glycerol backbones to have exactly one ester bond and exactly one phosphate connection.
    # In a genuine LPA:
    #   - One hydroxyl of glycerol is esterified with an acyl chain (identified by an OC(=O) substructure)
    #   - One other hydroxyl is directly bonded to the phosphate group.
    for glycerol in glycerol_matches:
        glycerol_atoms = set(glycerol)  # indices of atoms in the glycerol fragment

        # (a) Count ester bonds attached to the glycerol.
        ester_count = 0
        for o_idx in ester_oxygens:
            # If the oxygen that is part of the ester linkage is within the glycerol backbone, count it.
            if o_idx in glycerol_atoms:
                ester_count += 1

        if ester_count != 1:
            # This particular glycerol match is not monoacyl.
            continue

        # (b) Check that phosphate group is directly attached to this glycerol.
        # We expect that one of the phosphate oxygen atoms bonds to a carbon in the glycerol.
        connection_count = 0
        # Loop over each phosphate match.
        for phos_match in phos_matches:
            # In the phosphate SMARTS, the P atom is one of the indices.
            for idx in phos_match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 15:  # phosphorus atom
                    # Now look at its oxygen neighbors.
                    for nbr in atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 8:  # oxygen
                            # Check neighbors of the oxygen (ignoring the P atom) for a carbon that is in glycerol.
                            for nbr2 in nbr.GetNeighbors():
                                if nbr2.GetIdx() != atom.GetIdx() and nbr2.GetAtomicNum() == 6:
                                    if nbr2.GetIdx() in glycerol_atoms:
                                        connection_count += 1
                                        break  # found one connection from this O
                    # We only need one phosphate-to-glycerol connection.
                    if connection_count > 0:
                        break
            if connection_count > 0:
                break
        if connection_count != 1:
            # This glycerol match does not have a unique phosphate attachment;
            # it might be part of a larger headgroup.
            continue

        # (c) Check molecular weight.
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 250:
            return False, "Molecular weight too low to be an LPA."

        # If we have found a glycerol backbone with exactly one ester bond and exactly one phosphate connection,
        # then we consider the molecule a lysophosphatidic acid.
        return True, "Molecule contains a phosphate group directly linked to a glycerol backbone with exactly one acyl ester bond."
    
    # If none of the glycerol matches yield the proper pattern, then reject.
    return False, "No glycerol backbone found with exactly one acyl ester bond and one phosphate linkage."

# Example usage:
if __name__ == "__main__":
    # An example SMILES; you can test with any SMILES from the list.
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(O)(O)=O"  # PA(17:1(9Z)/0:0) expected to be LPA
    result, reason = is_lysophosphatidic_acids(test_smiles)
    print(result, reason)