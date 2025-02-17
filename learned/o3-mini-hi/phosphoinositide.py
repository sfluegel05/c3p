"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: Phosphoinositide
Definition: Any phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.
That is, the molecule must contain a glycerophosphatidylinositol “core” featuring two acyl chains (attached via ester bonds)
and a phosphate group that connects the glycerol backbone to an inositol ring where at least one hydroxyl is phosphorylated.
We now require that in addition to having ≥2 P atoms and ≥2 acyl ester (C(=O)O) groups,
a substructure representing a phosphorylated inositol head is present.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    
    Heuristic criteria:
      1. SMILES string must be valid.
      2. Molecule must have at least two phosphorus atoms.
      3. Molecule must contain at least two acyl ester groups (the C(=O)O motif).
      4. Molecule must contain a glycerophosphoinositol core. We require that a phosphate group
         (with zero formal charge) is connected via an oxygen to a cyclohexane ring bearing several hydroxyl groups
         (a common motif for inositol). The SMARTS pattern below looks for the fragment:
           O–P(=O)(O)O[C;R]1... 
         where the “1” ring is a cyclohexane ring having at least 4–5 hydroxyl substituents.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phosphoinositide, False otherwise.
        str: Reason explaining the classification decision.
    """
    
    # 1. Try to parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 2. Check for at least two phosphorus atoms (atomic number 15).
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 2:
        return False, f"Only {p_count} phosphorus atom(s) found; phosphoinositides require at least 2."
    
    # 3. Look for at least two acyl ester groups.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} acyl ester group(s); need at least 2 acyl chains."
    
    # 4. Look for a glycerophosphoinositol core.
    # First, define a SMARTS for a phosphorylated inositol head.
    # The phosphate fragment: O-P(=O)(O)O must be connected to a cyclohexane ring with several hydroxyl groups.
    # Here we combine two motifs:
    #    (a) The headgroup connectivity: OP(=O)(O)O[C;R]1 .... 1
    #    (b) The inositol ring: a cyclohexane ring (R6) where most ring carbons bear at least one -OH.
    #
    # This SMARTS looks for a phosphate attached via oxygen to a ring starting at atom index 1.
    # We require 4-6 oxygens on the ring; here, for simplicity, we demand at least 3 -OH groups.
    inositol_pattern = Chem.MolFromSmarts("OP(=O)(O)O[C;R]1[C;R]([OX2H])[C;R]([OX2H])[C;R]([OX2H])[C;R]([OX2H])[C;R]1[OX2H]")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        # If not found with the full pattern, try a looser pattern for an inositol ring:
        # Look for a cyclohexane ring bearing at least 3 hydroxyl groups.
        loose_inositol = Chem.MolFromSmarts("[$([C;R6]([OX2H])([OX2H])([OX2H]))]")
        ring_matches = mol.GetSubstructMatches(loose_inositol)
        if not ring_matches:
            return False, "No inositol ring (cyclohexane with multiple hydroxyls) attached to a phosphate detected."
        else:
            return False, "Inositol-like ring detected, but not attached via a neutral phosphate as required."
    else:
        # For each match, confirm that the phosphate involved is neutral.
        for match in inositol_matches:
            # In the inositol SMARTS, atom at position 1 is the oxygen attached to phosphorus,
            # atom 2 is the phosphorus. (Counting depends on how RDKit assigns indices; here we assume the
            # pattern: O (index0) - P (index1) ...).
            # We retrieve the phosphorus atom from the match.
            # Find the phosphorus among the atoms in the match.
            p_atom = None
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 15:
                    p_atom = atom
                    break
            if p_atom is not None and p_atom.GetFormalCharge() == 0:
                return True, ("Molecule contains at least two phosphorus atoms, at least two acyl ester groups, "
                              "and a glycerophosphoinositol core in which a neutral phosphate bridges a glycerol backbone "
                              "to an inositol ring.")
        return False, "Found a putative glycerophosphoinositol core, but the phosphate group is charged."
    
# Debug/test code (only executes when module is run directly)
if __name__ == "__main__":
    test_smiles = [
        # True positive example: PIP(18:0/16:0)
        "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O",
        # True positive example: 1-stearoyl-2-oleoyl-sn-glycero-3-phospho-1D-myo-inositol 4-phosphate
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC",
        # False positive candidate: similar motif but with anionic phosphate (should be rejected)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](OP([O-])([O-])=O)[C@H]1O)OC(=O)CCCCCCC\\C=C/CCCCCCCC"
    ]
    
    for smi in test_smiles:
        result, reason = is_phosphoinositide(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("----------")