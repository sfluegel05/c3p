"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: Phosphoinositide
Definition: Any phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.
That is, the molecule must contain a glycerophosphatidylinositol “core” (with two acyl chains attached via ester bonds
and a phosphate linking to an inositol ring) where at least one hydroxyl on the inositol has been phosphorylated.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.

    Heuristic criteria:
    1. The molecule must be a valid structure.
    2. It must contain at least two phosphorus atoms (one in the glycerol linkage and at least one extra on inositol).
    3. It must have at least two acyl ester groups (i.e. -C(=O)O- fragments) representing acyl chains.
    4. It must contain an inositol ring. Here we search for a six-membered ring of carbons where each ring atom 
       has at least one substituent oxygen (as would occur for hydroxyl/phosphate substituents) and where at least one of 
       these oxygens is phosphorylated (i.e. connected to a phosphorus).
       
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise.
        str: Explanation for the classification.
    """
    # 1. Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check the total number of phosphorus atoms – need at least 2.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 2:
        return False, f"Only {p_count} phosphorus atom(s) found; phosphoinositides need at least 2."
    
    # 3. Check for the presence of at least two acyl ester groups (i.e. a -C(=O)O- fragment).
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Only {len(ester_matches)} acyl ester group(s) found; need at least two as acyl chains."
    
    # 4. Look for an inositol ring:
    # Here we loop over all rings of size 6. For each such ring composed solely of carbons,
    # we check that every ring carbon has at least one substituent oxygen (not in the ring).
    # Furthermore, at least one substituent oxygen must be bound to a phosphorus.
    ring_info = mol.GetRingInfo().AtomRings()
    inositol_found = False
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Verify that all atoms in the ring are carbons.
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue

        all_have_oxygen_substituent = True
        phosphorylated_substituent = False
        # For each atom in the ring, check that it has at least one neighbor (not in the ring) that is oxygen.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            has_oxy = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                if neighbor.GetAtomicNum() == 8:  # oxygen
                    has_oxy = True
                    # Check if this oxygen is phosphorylated (i.e. connected to phosphorus)
                    for subneighbor in neighbor.GetNeighbors():
                        if subneighbor.GetAtomicNum() == 15:
                            phosphorylated_substituent = True
                            break
                    # No need to check further neighbors for this atom if one oxygen is found.
                    break
            if not has_oxy:
                all_have_oxygen_substituent = False
                break
        if all_have_oxygen_substituent and phosphorylated_substituent:
            inositol_found = True
            break

    if not inositol_found:
        return False, "No inositol ring detected that meets the substitution criteria (six-membered carbon ring with oxygen substituents and at least one phosphorylated substituent)."
    
    # If all three tests pass, we classify the molecule as a phosphoinositide.
    return True, "Molecule contains at least two phosphorus atoms, at least two acyl chains, and an inositol ring with phosphorylation."

# Debug/test code (if the module is run as a standalone program)
if __name__ == "__main__":
    test_smiles = [
        # Example: 1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-3,4-bisphosphate)
        "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCC",
        # Example: PIP(18:0/16:0)
        "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
    ]
    
    for smi in test_smiles:
        result, reason = is_phosphoinositide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")