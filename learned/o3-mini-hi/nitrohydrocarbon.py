"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C‐nitro compound that is a hydrocarbon in which one or more of the hydrogens
has been replaced by nitro groups.
Strategy:
  1. Parse the SMILES string and add explicit hydrogens.
  2. Find the nitro group using the SMARTS pattern [N+](=O)[O-].
  3. For each nitro group, check that at least one neighboring carbon (not part of another nitro group)
     shows evidence that a hydrogen has been replaced by the nitro group.
     • For aliphatic (sp3) carbons, compare the expected H count (4 - number of non-nitro neighbors) to the current H count.
     • For aromatic carbons, require that if the carbon is bonded only to aromatic neighbors (aside from the nitro),
       then its remaining H count is 0 rather than the usual 1.
  4. Remove all nitro group atoms from the molecule, and verify that every remaining heavy atom is either C or H.
  5. (Adjustment) We remove the requirement that a scaffold carbon must still have an explicit hydrogen.
     This allows fully substituted aromatic scaffolds to be classified as nitrohydrocarbons.
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
    # Parse SMILES and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define nitro group SMARTS pattern [N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"
    
    # Record indices of all atoms in nitro groups to later ignore them in scaffold check
    nitro_atom_indices = set()
    for match in nitro_matches:
        nitro_atom_indices.update(match)
    
    valid_substitution_found = False
    
    # For every nitro group, check that it is attached to at least one carbon that shows evidence
    # of having lost a hydrogen.
    for match in nitro_matches:
        # Identify the nitrogen atom (atomic number 7)
        nitro_n_idx = None
        for idx in match:
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7:
                nitro_n_idx = idx
                break
        if nitro_n_idx is None:
            continue  # Should not happen
        
        nitrogen = mol.GetAtomWithIdx(nitro_n_idx)
        carbon_attached = False
        for nbr in nitrogen.GetNeighbors():
            # Consider only carbon neighbors
            if nbr.GetAtomicNum() != 6:
                continue
            carbon_attached = True
            # Get neighbors of this carbon not belonging to any nitro group.
            non_nitro_neighbors = [n for n in nbr.GetNeighbors() if n.GetIdx() not in nitro_atom_indices]
            # Record the current hydrogen count on that carbon.
            current_H = nbr.GetTotalNumHs()
            if nbr.GetIsAromatic():
                # For an aromatic carbon on a benzene-like ring, typically a carbon is CH (i.e. 1 H)
                # If a nitro group replaces that one hydrogen then current_H should be 0.
                # Additionally, require that all non-nitro neighbors are aromatic.
                if all(n.GetIsAromatic() for n in non_nitro_neighbors):
                    # Note: In polyaromatics, the normal expected hydrogen count may be 0 if fully substituted,
                    # so we only flag as valid substitution when we see the specific pattern of substitution.
                    # Here we require exactly two aromatic neighbors and 0 hydrogens.
                    if len(non_nitro_neighbors) == 2 and current_H == 0:
                        valid_substitution_found = True
            else:
                # For an aliphatic carbon, normally the ideal number of hydrogens would be (4 - number of heavy neighbors).
                # Since one H is replaced by the nitro group, we expect current_H == (ideal - 1).
                d = len(non_nitro_neighbors)
                ideal_H = 4 - d
                if ideal_H > 0 and current_H == ideal_H - 1:
                    valid_substitution_found = True
        if not carbon_attached:
            return False, "A nitro group is not attached to any carbon atom"
    
    if not valid_substitution_found:
        return False, ("No nitro group appears to be replacing a hydrogen on a carbon. "
                       "The attachment carbons might be too substituted for simple hydrogen replacement.")
    
    # Remove nitro-group atoms and check that the remaining scaffold is a pure hydrocarbon (only C and H)
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_indices:
            continue
        if atom.GetAtomicNum() not in (1, 6):
            return False, f"Scaffold contains non-hydrocarbon atom: {atom.GetSymbol()}"
    
    return True, "Molecule is a nitrohydrocarbon: It is a hydrocarbon with at least one nitro group replacing a hydrogen."

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
        "4-[(1-hydroxy-2-naphthalenyl)hydrazinylidene]-7-nitro-3-oxo-1-naphthalenesulfonic acid": "C1=CC=C2C(=C1)C=CC(=C2O)NN=C3C4=C(C=C(C=C4)[N+](=O)[O-])C(=CC3=O)S(=O)(=O)O",
        "3-(5,8-dinitro-1,3-dioxo-2-benzo[de]isoquinolinyl)benzoic acid": "C1=CC(=CC(=C1)N2C(=O)C3=CC(=CC4=CC(=CC(=C43)C2=O)[N+](=O)[O-])[N+](=O)[O-])C(=O)O",
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