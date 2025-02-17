"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: Phosphoinositide
Definition: Any phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.
That is, the molecule must contain a glycerophosphatidylinositol “core” featuring two acyl chains (attached via ester bonds)
and a phosphate group that connects the glycerol backbone to an inositol ring where at least one of the inositol hydroxyls is phosphorylated.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.

    Heuristic criteria:
      1. Molecule parses correctly.
      2. Contains at least two phosphorus atoms.
      3. Has at least two acyl ester groups (the C(=O)O motif) indicating acyl chains.
      4. Contains a glycerophosphoinositol core. For this we look for a specific connectivity:
         an –O–C–C–O–P(=O)(O)O fragment that connects a glycerol-like moiety to the inositol head.
         Also we require that the phosphorus in this fragment is neutral.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phosphoinositide, False otherwise.
        str: Explanation for the classification.
    """
    # 1. Parse SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 2. Check total phosphorus count (atomic number 15) should be at least 2.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 2:
        return False, f"Found only {p_count} phosphorus atom(s); phosphoinositides require at least 2."
    
    # 3. Check for at least two acyl ester groups.
    # The ester (acyl chain) pattern is defined as a C(=O)O fragment.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} acyl ester group(s); need at least 2 acyl chains."
    
    # 4. Look for a glycerophosphoinositol core.
    # We define a substructure pattern for the key connectivity:
    # An oxygen connected to two carbons then an oxygen connected to a phosphorus with three oxygens.
    # The SMARTS pattern below is written to capture an –O–C–C–O–P(=O)(O)O fragment.
    # Note: Stereochemistry or charge markers are left out so that we can match most common variants.
    gp_pattern = Chem.MolFromSmarts("[O]CCOP(=O)(O)O")
    gp_matches = mol.GetSubstructMatches(gp_pattern)
    if not gp_matches:
        return False, "Glycerophosphoinositol core (O–C–C–O–P(=O)(O)O fragment) not detected."
    
    # Check that in at least one match the phosphate group is neutral
    valid_core_found = False
    for match in gp_matches:
        # In the pattern [O]CCOP(=O)(O)O, the phosphorus should be at index 4.
        p_atom = mol.GetAtomWithIdx(match[4])
        if p_atom.GetFormalCharge() == 0:
            valid_core_found = True
            break
    if not valid_core_found:
        return False, "Glycerophosphoinositol core detected but the phosphate group is not neutral."
    
    # If all tests pass, we consider it a phosphoinositide.
    return True, "Molecule contains at least two phosphorus atoms, at least two acyl ester groups, and a neutral glycerophosphoinositol core linking to an inositol ring."

# Debug/test code (when run standalone)
if __name__ == "__main__":
    test_smiles = [
        # Example: PIP(18:0/16:0)
        "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O",
        # Example: 1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol 4-phosphate
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC",
        # Example: a false positive candidate (with anionic phosphate in the glycerol core)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP([O-])([O-])=O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC"
    ]
    
    for smi in test_smiles:
        result, reason = is_phosphoinositide(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("----------")