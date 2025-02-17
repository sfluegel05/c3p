"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate

Definition:
  An sn-glycero-3-phosphate compound having a single unspecified acyl group
  at either position 1 or position 2. That means that the glycerol backbone 
  (with three carbons, where carbon-3 is phosphorylated) must have exactly one 
  acyl (ester) substitution on either sn-1 or sn-2 and the remaining hydroxyl should be free.
  
Approach:
  1. Parse the SMILES string.
  2. Require exactly one phosphorus atom and reject molecules containing nitrogen.
  3. Add explicit hydrogens (so that free –OH groups become visible).
  4. Look for a glycerol-3-phosphate backbone by finding a substructure pattern 
     where three (mapped) carbons connect linearly with the third attached to phosphate.
     (We use the pattern: "[C:1][C:2][C:3]OP(=O)(O)O" where mapped atoms 1,2,3 denote sn-1, sn-2, sn-3.)
  5. For the two “glycerol” carbons that are not connected to phosphate, inspect the 
     oxygen substituents. One must be acylated (ester linkage to a carbonyl group) and the other must be a free hydroxyl (only bound to hydrogen apart from the glycerol carbon).
  6. To decide if an oxygen is acylated we require that it is connected to a carbon (the carbonyl carbon)
     that in turn is double-bonded to an oxygen.
  7. Finally, if exactly one of sn-1 and sn-2 is acylated, report that the molecule is of the correct class.
  
If any test fails we return (False, <reason>).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True if the molecule is classified as a monoacyl-sn-glycerol 3-phosphate,
                     False otherwise, together with an explanation.
    """
    # 1. Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Require exactly one phosphorus atom
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, f"Molecule must contain exactly one phosphorus atom (found {len(phosphorus_atoms)})"
    
    # Reject molecules containing nitrogen
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if nitrogen_atoms:
        return False, "Molecule contains nitrogen atoms, indicating an alternative headgroup"
    
    # 3. Add explicit hydrogens so that we can detect free hydroxyl groups.
    mol = Chem.AddHs(mol)
    
    # 4. Look for a glycerol-3-phosphate backbone.
    # We expect a linear three-carbon fragment where C3 has an O-bound phosphate.
    # Note: our mapped pattern only returns indices for the mapped atoms.
    backbone_smarts = "[C:1][C:2][C:3]OP(=O)(O)O"
    backbone_pattern = Chem.MolFromSmarts(backbone_smarts)
    if backbone_pattern is None:
        return False, "Invalid SMARTS for glycerol-phosphate backbone"
    
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if not backbone_matches:
        return False, "No glycerol-phosphate backbone found"
    
    # For simplicity, use the first backbone match.
    # According to the pattern, match[0] = sn-1 carbon, match[1] = sn-2 carbon, match[2] = sn-3 carbon.
    sn1_idx, sn2_idx, sn3_idx = backbone_matches[0]
    sn1_atom = mol.GetAtomWithIdx(sn1_idx)
    sn2_atom = mol.GetAtomWithIdx(sn2_idx)
    # We do not need to inspect sn3 beyond the attachment to phosphate.
    
    # Helper function: given an oxygen atom connected to a glycerol carbon, decide if it is acylated.
    # We consider an oxygen "acylated" if it is bound (via a single bond) to a carbon that
    # has at least one double bond to an oxygen.
    def is_acylated(oxy_atom):
        # Find the neighbor (other than the glycerol carbon) that is carbon.
        for neigh in oxy_atom.GetNeighbors():
            # Look for a carbon neighbor not part of the glycerol backbone.
            if neigh.GetAtomicNum() == 6:
                # Check bonds of the carbon for a C=O
                for nb in neigh.GetNeighbors():
                    if nb.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(neigh.GetIdx(), nb.GetIdx())
                        if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            return True
        return False

    # Check substituents on sn1 and sn2 that are oxygens not involved in the backbone link.
    # For each glycerol carbon, ignore bonds to carbons that are part of the backbone.
    def get_external_oxygen(carbon_atom, backbone_neighbors):
        # Return a list of oxygen neighbors that are not in backbone_neighbors.
        oxygens = []
        for neigh in carbon_atom.GetNeighbors():
            if neigh.GetAtomicNum() == 8 and neigh.GetIdx() not in backbone_neighbors:
                oxygens.append(neigh)
        return oxygens

    # For sn1, backbone neighbors are:
    sn1_backbone = {sn2_idx}
    # For sn2, backbone neighbors:
    sn2_backbone = {sn1_idx, sn3_idx}
    
    sn1_oxygens = get_external_oxygen(sn1_atom, sn1_backbone)
    sn2_oxygens = get_external_oxygen(sn2_atom, sn2_backbone)
    
    # Count how many of these oxygens are acylated.
    acyl_count = 0
    acyl_positions = []  # keep track of which position is acylated
    # Also check that the free OH has an explicit hydrogen.
    free_oh_ok = [False, False]  # for sn1 and sn2 respectively
    
    # For sn1:
    sn1_acyl = False
    sn1_free = False
    for oxy in sn1_oxygens:
        if is_acylated(oxy):
            sn1_acyl = True
        else:
            # Check that this oxygen is a free hydroxyl. It must have at least one H.
            # Counting explicit hydrogens.
            n_H = oxy.GetTotalNumHs(includeNeighbors=True)
            if n_H > 0:
                sn1_free = True
    # For sn2:
    sn2_acyl = False
    sn2_free = False
    for oxy in sn2_oxygens:
        if is_acylated(oxy):
            sn2_acyl = True
        else:
            n_H = oxy.GetTotalNumHs(includeNeighbors=True)
            if n_H > 0:
                sn2_free = True

    # We now require that exactly one of the two positions is acylated and the other presents a free OH.
    if sn1_acyl and sn2_acyl:
        return False, "Both sn-1 and sn-2 appear acylated; only one acyl group is allowed"
    if (not sn1_acyl) and (not sn2_acyl):
        return False, "Neither sn-1 nor sn-2 is acylated; exactly one acyl group must be present"
    if sn1_acyl and (not sn2_free):
        return False, "sn-2 oxygen does not appear as a free hydroxyl (expected a free OH group)"
    if sn2_acyl and (not sn1_free):
        return False, "sn-1 oxygen does not appear as a free hydroxyl (expected a free OH group)"
    
    # If we reached this point, we have exactly one acyl substitution at either sn-1 or sn-2.
    return True, "Molecule has a glycerol-phosphate backbone with exactly one acyl ester group at sn-1 or sn-2"

    
# For testing purposes, below are some candidate SMILES from the provided list.
if __name__ == "__main__":
    test_smiles_list = [
        "CCCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O",   # 1-icosanoyl-sn-glycero-3-phosphate (expected True)
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCC)(O)(O)=O",          # 1-nonadecanoyl-glycero-3-phosphate (expected True)
        "P(OC[C@H](O)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)(O)(O)=O",  # 1-(6Z,9Z,12Z,15Z-octadecatetraenoyl)-... (expected True)
        "O(C[C@H](COP(=O)(O)O)O)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC",   # 1-linoleoyl-sn-glycerol 3-phosphate (expected True)
        "CCCCCCCCCCCCCCCOC[C@H](COP(O)(O)=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC",  # 1-palmityl-2-arachidonoyl-sn-glycero-3-phosphate (expected False)
    ]
    for smi in test_smiles_list:
        res, msg = is_monoacyl_sn_glycerol_3_phosphate(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {msg}\n")