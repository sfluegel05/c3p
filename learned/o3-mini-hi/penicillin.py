"""
Classifies: CHEBI:17334 penicillin
"""
#!/usr/bin/env python
"""
Classifies: Penicillins
Definition:
  “Any member of the group of substituted penams containing two methyl substituents at position 2,
   a carboxylate substituent at position 3 and a carboxamido group at position 6.”
   
The detection is carried out in parts:
  1) Identify a fused bicyclic penam core by:
       - Finding a four‐membered ring (β‐lactam) with a nitrogen and with at least one carbonyl group.
       - Finding a five‐membered ring that shares exactly two atoms with the 4-membered ring and 
         contains a sulfur (thiazolidine).
  2) Verify that the core bears two methyl substituents at the position adjacent to sulfur.
  3) Check that a carboxylate substituent is present (either C(=O)O or C(=O)[O-]).
  4) Check that a carboxamido fragment (N–C(=O)) is present.
  
Note:
  This implementation uses both ring‐analysis using RDKit’s ring info and simple SMARTS matches.
  Because penicillin related molecules show variable substituents and stereochemistry, we try to be flexible.
"""

from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a penicillin.
    A molecule is considered a penicillin if:
      - It contains a fused bicyclic penam core. That is, a four‐membered β‐lactam ring fused with a five‐membered thiazolidine;
      - The thiazolidine ring bears two methyl groups on the carbon adjacent to the sulfur (position 2);
      - It contains at least one carboxylate substituent (C(=O)O or C(=O)[O-]);
      - It contains at least one carboxamido fragment (N–C(=O), for example at position 6).
      
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule appears to be a penicillin, otherwise False.
        str: Explanation of the decision.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Identify the penicillin fused bicyclic core ---
    # We define a helper function that returns True if the molecule has a fused
    # bicyclic penam core (a four-membered beta-lactam ring fused with a five-membered thiazolidine)
    def has_penicillin_core(mol):
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()  # list of tuples of atom indices
        # Look for a 4-membered ring that has a nitrogen and a carbonyl moiety
        for ring4 in rings:
            if len(ring4) != 4:
                continue
            ring4_set = set(ring4)
            atoms_ring4 = [mol.GetAtomWithIdx(idx) for idx in ring4_set]
            has_nitrogen = any(atom.GetSymbol() == "N" for atom in atoms_ring4)
            has_carbonyl = False
            for atom in atoms_ring4:
                if atom.GetSymbol() == "C":
                    for bond in atom.GetBonds():
                        # check if bond is double and the other atom is oxygen
                        if bond.GetBondType().name == "DOUBLE":
                            other = bond.GetOtherAtom(atom)
                            if other.GetSymbol() == "O":
                                has_carbonyl = True
                                break
                    if has_carbonyl:
                        break
            if not (has_nitrogen and has_carbonyl):
                continue
            # Now check for a 5-membered ring that fuses with the found 4-membered ring:
            for ring5 in rings:
                if len(ring5) != 5:
                    continue
                common = set(ring5).intersection(ring4_set)
                if len(common) == 2:  # fused ring sharing two atoms
                    # Check that the 5-membered ring contains a sulfur atom
                    atoms_ring5 = [mol.GetAtomWithIdx(idx) for idx in ring5]
                    if any(atom.GetSymbol() == "S" for atom in atoms_ring5):
                        return True
        return False

    if not has_penicillin_core(mol):
        return False, "Molecule does not contain the expected fused penam core (β-lactam fused to thiazolidine)"
    
    # --- Step 2: Check for two methyl substituents at the key position (position 2) ---
    # We look for an S-bound sp3 carbon bearing two CH3 groups.
    # The SMARTS below is written to detect an sp3 carbon (C) directly bonded to an S, that in turn carries two CH3 groups.
    dimethyl_pattern = Chem.MolFromSmarts("[SX2]-[CX4]([CH3])([CH3])")
    if dimethyl_pattern is None:
        return False, "Error in dimethyl SMARTS"
    if not mol.HasSubstructMatch(dimethyl_pattern):
        return False, "Missing the two methyl substituents at the C adjacent to the sulfur (position 2)"
    
    # --- Step 3: Check for the carboxylate substituent ---
    # We accept either deprotonated or acid forms.
    carboxylate_smarts1 = Chem.MolFromSmarts("C(=O)[O-]")
    carboxylate_smarts2 = Chem.MolFromSmarts("C(=O)O")
    if carboxylate_smarts1 is None or carboxylate_smarts2 is None:
        return False, "Error in carboxylate SMARTS"
    if not (mol.HasSubstructMatch(carboxylate_smarts1) or mol.HasSubstructMatch(carboxylate_smarts2)):
        return False, "Missing a carboxylate substituent (C(=O)O or C(=O)[O-])"
    
    # --- Step 4: Check for the carboxamido group ---
    # We look for a fragment: N-C(=O). (This is a simplified check.)
    carboxamido_smarts = Chem.MolFromSmarts("N-C(=O)")
    if carboxamido_smarts is None:
        return False, "Error in carboxamido SMARTS"
    if not mol.HasSubstructMatch(carboxamido_smarts):
        return False, "Missing a carboxamido group (N-C(=O) fragment)"
    
    return True, "Molecule has a penicillin fused core plus the required substituents (dimethyl at pos.2, carboxylate at pos.3, carboxamido at pos.6)"

# Example usage:
if __name__ == "__main__":
    # Example: Penicillin K SMILES (one of many examples)
    test_smiles = "CCCCCCCC(=O)N[C@H]1[C@H]2SC(C)(C)[C@@H](N2C1=O)C(O)=O"
    result, reason = is_penicillin(test_smiles)
    print("Is penicillin?", result)
    print("Reason:", reason)