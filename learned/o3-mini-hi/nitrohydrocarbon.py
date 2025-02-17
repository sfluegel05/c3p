"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C‐nitro compound that is a hydrocarbon in which one or more of the hydrogens
has been replaced by nitro groups.

Our strategy is as follows:
  1. Parse the SMILES and add explicit hydrogens.
  2. Look for the nitro group using a SMARTS pattern: [N+](=O)[O-].
  3. For each nitro group, ensure that its nitrogen is bonded to at least one carbon.
  4. For at least one nitro group, require evidence that the carbon (if the nitro were replaced
     by hydrogen) would gain an extra hydrogen relative to its current state.
       • For an aliphatic carbon: In a saturated hydrocarbon, a carbon with d heavy neighbors
         normally would have (4 - d) hydrogens. If a nitro replaces a hydrogen, the current H count
         should be (4 - d - 1) = 3 - d.
       • For an aromatic carbon: In benzene a C normally is CH (i.e. 1 H) when bonded to exactly
         two other aromatic carbons. If the nitro replaces that hydrogen then the current H count
         should be 0. (We will also require that every non‐nitro neighbor of the aromatic C is itself
         aromatic so that substitution really took place on the uniform benzene scaffold.)
  5. Finally, remove (ignore) all nitro group atoms and verify that every remaining heavy atom is
     carbon.
  6. Also check that at least one scaffold carbon (i.e. not part of any nitro group) still carries
     some hydrogen(s) (so that the scaffold isn’t fully substituted by non‐hydrogen groups).

If any of these conditions fails then the molecule is not classified as a nitrohydrocarbon.
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule matches the nitrohydrocarbon criteria, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES and add hydrogens (so that we can examine explicit H counts)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define nitro group SMARTS pattern [N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"
    
    # Record indices of atoms that belong to a nitro group.
    nitro_atom_indices = set()
    for match in nitro_matches:
        nitro_atom_indices.update(match)
    
    # Flag to see if we have found at least one nitro substitution on a carbon that “lost” a hydrogen.
    valid_substitution_found = False
    
    # For each nitro group match...
    for match in nitro_matches:
        # Identify the nitrogen atom in this nitro group (atomic number 7)
        nitro_n_idx = None
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                nitro_n_idx = idx
                break
        if nitro_n_idx is None:
            continue  # should not happen
        
        nitrogen = mol.GetAtomWithIdx(nitro_n_idx)
        # Look for at least one neighboring carbon (which is not part of a nitro group)
        carbon_found = False
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            carbon_found = True
            # Examine this attachment carbon.
            # Get its neighbors that are not part of any nitro group.
            non_nitro_neighbors = [n for n in nbr.GetNeighbors() if n.GetIdx() not in nitro_atom_indices]
            # Note the current explicit hydrogen count.
            current_H = nbr.GetTotalNumHs()
            if nbr.GetIsAromatic():
                # For an aromatic carbon (we require that all non-nitro neighbors are aromatic)
                # In benzene a carbon (CH) normally is bonded to two other aromatic carbons.
                if any(not x.GetIsAromatic() for x in non_nitro_neighbors):
                    # If any non-nitro neighbor is not aromatic, assume it is already substituted.
                    continue
                # For a benzene carbon, if no substitution occurred it should have 1 hydrogen.
                # If a nitro group has replaced that H, then current_H should be 0.
                if len(non_nitro_neighbors) == 2 and current_H == 0:
                    valid_substitution_found = True
            else:
                # For an aliphatic (sp3) carbon:
                # In a fully saturated hydrocarbon, an unsubstituted carbon would have: 
                #   ideal_H = 4 - (number of heavy (non-H) neighbors)
                # When a nitro group replaces a hydrogen, current_H = ideal_H - 1.
                d = len(non_nitro_neighbors)  # current heavy (non-nitro) substituents
                ideal_H = 4 - d
                if ideal_H > 0 and current_H == ideal_H - 1:
                    valid_substitution_found = True
        if not carbon_found:
            return False, "A nitro group is not attached to any carbon atom"
    
    if not valid_substitution_found:
        return False, ("No nitro group appears to be replacing a hydrogen on a carbon. "
                       "The attachment carbons are too substituted to be simple hydrogen replacements.")
    
    # Check that when we ignore all nitro group atoms the scaffold is a pure hydrocarbon (only C and H).
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_indices:
            continue
        # Allow carbon (atomic number 6) and hydrogen (atomic number 1) only.
        if atom.GetAtomicNum() not in (1, 6):
            return False, f"Scaffold contains non-hydrocarbon atom: {atom.GetSymbol()}"
    
    # Additionally, ensure that at least one scaffold carbon (not in a nitro group) still carries hydrogen.
    scaffold_has_H = False
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_indices:
            continue
        if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() > 0:
            scaffold_has_H = True
            break
    if not scaffold_has_H:
        return False, "No hydrogen remains on the hydrocarbon scaffold"
    
    return True, "Molecule is a nitrohydrocarbon: nitro groups are attached in place of hydrogen on a pure hydrocarbon scaffold"

# Example usage / testing:
if __name__ == "__main__":
    test_examples = {
        "3,7-Dinitrofluoranthene": "[O-][N+](=O)c1cccc-2c1-c1cccc3c(ccc-2c13)[N+]([O-])=O",
        "1-nitroheptane": "[O-][N+](=O)CCCCCCC",
        "3,9-Dinitrofluoranthene": "[O-][N+](=O)c1ccc-2c(c1)-c1ccc([N+]([O-])=O)c3cccc-2c13",
        "1-nitropyrene": "[O-][N+](=O)c1ccc2ccc3cccc4ccc1c2c34",
        "nitroethene": "[O-][N+](=O)C=C",
        "2,7-dinitronaphthalene": "[O-][N+](=O)c1ccc2ccc(cc2c1)[N+]([O-])=O",
        "3-Nitrofluoranthene": "[O-][N+](=O)c1ccc2-c3ccccc3-c3cccc1c23",
        "1-nitropropane": "CCC[N+]([O-])=O",
        "2-nitronaphthalene": "[O-][N+](=O)c1ccc2ccccc2c1",
        "1,5-dinitronaphthalene": "[O-][N+](=O)c1cccc2c(cccc12)[N+]([O-])=O",
        "1,4-dinitronaphthalene": "[O-][N+](=O)c1ccc([N+]([O-])=O)c2ccccc12",
        "2-nitrofluorene": "[O-][N+](=O)c1ccc-2c(Cc3ccccc-23)c1",
        "nitrobenzene": "[O-][N+](=O)c1ccccc1",
        "1,3-dinitronaphthalene": "[O-][N+](=O)c1cc([N+]([O-])=O)c2ccccc2c1",
        # Some expected false positives:
        "6-Nitrobenzo[a]pyrene (expected False)": "[O-][N+](=O)c1c2ccccc2c2ccc3cccc4ccc1c2c34",
        "4-nitrotoluene (expected False)": "Cc1ccc(cc1)[N+]([O-])=O",
        "4-Nitrobiphenyl (expected False)": "[O-][N+](=O)c1ccc(cc1)-c1ccccc1",
        "6-Nitrochrysene (expected False)": "[O-][N+](=O)c1cc2c3ccccc3ccc2c2ccccc12",
        "cis-6-Nitro-p-mentha-1(7),2-diene (expected False)": "[O-][N+](=O)C1CC(C(C)C)C=CC1=C",
        "17beta-Nitro-5alpha-androstane (expected False)": "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCCC[C@]34C)[C@@H]1CC[C@@H]2[N+]([O-])=O",
        "2,4,6-trinitrotoluene (expected False)": "Cc1c(cc(cc1[N+]([O-])=O)[N+]([O-])=O)[N+]([O-])=O",
        "3,4-dinitrotoluene (expected False)": "C1=C(C(=CC=C1C)[N+]([O-])=O)[N+](=O)[O-]",
        "3-nitrotoluene (expected False)": "Cc1cccc(c1)[N+]([O-])=O",
        "2-methyl-1,3,5-trinitrocyclohex-2-en-1-ide (expected False)": "C=1(C)[C-]([N+](=O)[O-])CC(CC1[N+]([O-])=O)[N+](=O)[O-]",
        "1-nitrobutane (expected False)": "[O-][N+](=O)CCCC"
    }
    for name, smi in test_examples.items():
        result, reason = is_nitrohydrocarbon(smi)
        print(f"{name}: {result}. Reason: {reason}")