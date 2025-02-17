"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C‐nitro compound that is a hydrocarbon in which one or more of the hydrogens 
has been replaced by nitro groups.
In our approach the molecule is classified as a nitrohydrocarbon if:
  1. It contains one or more nitro groups (SMARTS "[N+](=O)[O-]"),
  2. Each nitro group is attached to at least one carbon atom,
  3. At least one of those attachment carbons would have originally had a hydrogen 
     (i.e. when the nitro group is “replaced” by hydrogen the heavy‐atom count on that carbon 
     is lower than the maximal allowed for that hybridization: 3 for aromatic C or 4 for aliphatic C),
  4. And when nitro groups are removed, every remaining heavy atom in the scaffold is carbon 
     and at least one carbon in that scaffold has at least one hydrogen.
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    Our method does the following:
       (1) Verify that at least one [N+](=O)[O-] group is present.
       (2) For each nitro group, ensure its nitrogen is bonded to at least one carbon.
       (3) Confirm that for at least one nitro group, the carbon to which it is attached
           shows evidence of “losing a hydrogen” (i.e. it is not fully substituted by heavy atoms).
       (4) Verify that if you ignore all the atoms that belong to nitro groups, every remaining 
           heavy atom (non‐hydrogen) is carbon.
       (5) Check that at least one scaffold carbon (i.e. not belonging to a nitro group) carries hydrogen(s).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule fulfills our criteria of being a nitrohydrocarbon, False otherwise.
        str: Explanation of the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find nitro groups using the SMARTS pattern for nitro: [N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"

    # Record all atoms that belong to any nitro group (indices)
    nitro_atom_indices = set()
    for match in nitro_matches:
        nitro_atom_indices.update(match)
    
    # For each nitro group, check that the N atom is attached to at least one carbon.
    # Also, for at least one nitro we want to see that the attached carbon was originally a CH.
    at_least_one_valid_substitution = False
    for match in nitro_matches:
        # Identify the nitrogen atom in the match (atomic number 7)
        nitro_n_idx = None
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                nitro_n_idx = idx
                break
        if nitro_n_idx is None:
            continue  # should not happen
        nitrogen = mol.GetAtomWithIdx(nitro_n_idx)
        # Look among its neighbors for a carbon atom.
        carbon_found = False
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_found = True
                # Now check if this carbon originally had a hydrogen available.
                # We do this by counting how many heavy (non-H) neighbors it has,
                # excluding any neighbor that is part of a nitro group.
                # For an aromatic carbon, ideally it should have 3 bonds (2 to carbons in the ring and 1 to H)
                # and for an aliphatic (sp3), ideally 4 bonds.
                nbr_indices = [n.GetIdx() for n in nbr.GetNeighbors() 
                               if n.GetAtomicNum() > 1 and n.GetIdx() not in nitro_atom_indices]
                if nbr.GetIsAromatic():
                    max_allowed = 3
                else:
                    max_allowed = 4
                # If the number of heavy neighbors (excluding nitro atoms) is less than the maximum,
                # then that carbon would have had at least one hydrogen.
                if len(nbr_indices) < max_allowed:
                    at_least_one_valid_substitution = True
        if not carbon_found:
            return False, "A nitro group is not attached to any carbon atom"
    
    if not at_least_one_valid_substitution:
        return False, "No nitro group appears to replace a hydrogen on a carbon (attachment carbon is fully substituted)"
    
    # Check the hydrocarbon scaffold: every heavy atom not in a nitro group must be a carbon.
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_indices:
            continue
        if atom.GetAtomicNum() != 6:
            return False, f"Scaffold contains non-hydrocarbon atom: {atom.GetSymbol()}"
    
    # Additionally, ensure that at least one scaffold carbon has a hydrogen.
    scaffold_has_H = False
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_indices:
            continue
        if atom.GetTotalNumHs() > 0:
            scaffold_has_H = True
            break
    if not scaffold_has_H:
        return False, "No hydrogen remains on the hydrocarbon scaffold"
    
    return True, "Molecule is a nitrohydrocarbon: nitro groups are attached in place of hydrogen on a pure hydrocarbon scaffold"

# Example usage / testing
if __name__ == "__main__":
    test_examples = {
        "3,7-Dinitrofluoranthene": "[O-][N+](=O)c1cccc-2c1-c1cccc3c(ccc-2c13)[N+]([O-])=O",
        "1-nitroheptane": "[O-][N+](=O)CCCCCCC",
        "Nitrobenzene": "[O-][N+](=O)c1ccccc1",
        "3,9-Dinitrofluoranthene": "[O-][N+](=O)c1ccc-2c(c1)-c1ccc([N+]([O-])=O)c3cccc-2c13",
        "1-nitropropane": "CCC[N+]([O-])=O",
        "4-nitrotoluene (expected False)": "Cc1ccc(cc1)[N+]([O-])=O",  # should be rejected if the substituted carbon was already substituted by CH3
        "6-Nitrobenzo[a]pyrene (expected False)": "[O-][N+](=O)c1c2ccccc2c2ccc3cccc4ccc1c2c34"
    }
    
    for name, smi in test_examples.items():
        result, reason = is_nitrohydrocarbon(smi)
        print(f"{name}: {result}. Reason: {reason}")