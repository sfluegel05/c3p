"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C‐nitro compound that is a hydrocarbon in which one or more of the hydrogens
has been replaced by nitro groups.
Our strategy (heuristic):
  1. Parse the SMILES and add explicit hydrogens.
  2. Look for nitro groups via the SMARTS pattern "[N+](=O)[O-]".
  3. For each nitro group, find the carbon atom that it is attached to.
     For that carbon we “simulate” what its hydrogen count would be if the –NO₂ group were replaced by H.
     • For an sp³ carbon, we assume that in the parent hydrocarbon the ideal hydrogen count is: 4 – (# non-nitro heavy neighbors).
     • For an sp² (non‐aromatic) carbon, we assume an ideal count of 3 – (# non-nitro heavy neighbors).
     • For an aromatic carbon (we expect a benzene‐like ring), a carbon that normally bears an H should have exactly one H;
       so if it gets substituted then its current explicit hydrogen count is 0.
     In our “simulation” we add back (virtually) one hydrogen per nitro group attached.
  4. Require that at least one nitro group is “valid” in that the effective H count of the attached C equals the expected number.
  5. Then “remove” all nitro-group atoms and verify that every remaining heavy atom is C (i.e. the remaining scaffold is a hydrocarbon).
If all these checks pass we classify the molecule as a nitrohydrocarbon.
Note: There is considerable ambiguity when a nitro group is attached to a carbon that is already substituted.
In these borderline cases our algorithm may err. This is one plausible improvement over a previous attempt.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon (a C‐nitro hydrocarbon) based on its SMILES.
    
    A nitrohydrocarbon is defined here as a molecule containing at least one nitro group ([N+](=O)[O-])
    that appears to have replaced one hydrogen on a carbon in the parent hydrocarbon.
    In addition, when the nitro-group atoms are removed the remaining scaffold must consist only of
    carbon and hydrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is judged to be a nitrohydrocarbon, False otherwise.
        str: Explanation of the decision.
    """
    # Parse SMILES and add explicit Hs.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define the nitro group pattern.
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"
    
    # Record all atom indices that belong to any nitro group.
    nitro_atom_indices = set()  # will contain indices of N and its two O's in each NO2
    # Also record a mapping from an attachment carbon index to the count of nitro groups attached.
    attachment_dict = {}  # key: carbon atom idx; value: number of nitro groups attached
    
    for match in nitro_matches:
        # Each match is a tuple of indices corresponding to the [N+](=O)[O-] pattern.
        # We now want to know which carbon(s) these nitro groups are attached to.
        # In a nitro group, only the nitrogen atom (atomic number 7) can bond to a carbon.
        nitro_n_idx = None
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                nitro_n_idx = idx
                break
        if nitro_n_idx is None:
            continue  # safety check
        # Mark all nitro atoms so they are ignored when assessing the scaffold.
        nitro_atom_indices.update(match)
        # Look at neighbors of the nitro nitrogen. We require at least one to be carbon.
        found_carbon = False
        nitrogen = mol.GetAtomWithIdx(nitro_n_idx)
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                found_carbon = True
                cidx = nbr.GetIdx()
                # Count the nitro group attachment for each carbon.
                attachment_dict[cidx] = attachment_dict.get(cidx, 0) + 1
        if not found_carbon:
            return False, "A nitro group is not attached to any carbon atom"
    
    # Check that at least one nitro group occurs on a carbon in a manner consistent with hydrogen replacement.
    valid_substitution_found = False
    for cidx, n_nitro in attachment_dict.items():
        carbon = mol.GetAtomWithIdx(cidx)
        # Count neighbors that are not part of any nitro group.
        non_nitro_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetIdx() not in nitro_atom_indices]
        d = len(non_nitro_neighbors)
        current_H = carbon.GetTotalNumHs()  # number of explicit H attached
        effective_H = current_H + n_nitro  # if we “restore” one H for each nitro group attached
        
        # Determine expected number of hydrogens on the carbon if it were in the original hydrocarbon.
        # For sp3 (non-aromatic) carbons we assume an ideal total of 4 bonds;
        # for non‐aromatic sp2, ideal total of 3; for aromatic carbons we assume a benzene type carbon:
        #   if it normally bears H, then expected H == 1 (if only two heavy neighbors) and if fully substituted (3 heavy neighbors) then H==0.
        if carbon.GetIsAromatic():
            # For simplicity, if the carbon has exactly 2 neighbors (in the de-nitro scaffold) we expect 1 H.
            # If it has 3 or more, we expect it to be fully substituted (0 H). (This rule works for benzene-like rings.)
            expected_H = 1 if d == 2 else 0
        else:
            # Determine ideal total bonds based on hybridization:
            hyb = carbon.GetHybridization()
            if hyb == rdchem.HybridizationType.SP2:
                ideal = 3
            else:
                ideal = 4
            expected_H = ideal - d
        
        # For debugging, one might print the values:
        # print(f"Carbon idx {cidx}: current_H={current_H}, attached_nitro={n_nitro}, effective_H={effective_H}, non-nitro neighbors={d}, expected_H={expected_H}")
        
        if effective_H == expected_H:
            valid_substitution_found = True
            # Once a valid substitution is found, we can break out.
            break

    if not valid_substitution_found:
        return False, ("No nitro group appears to have replaced an available hydrogen on a carbon. "
                       "The attachment site(s) may be too substituted for a simple H replacement.")
    
    # Remove (ignore) atoms belonging to nitro groups and ensure that every remaining heavy atom is C.
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_indices:
            continue
        if atom.GetAtomicNum() not in (1, 6):  # allow only hydrogen or carbon
            return False, f"Scaffold contains non-hydrocarbon atom: {atom.GetSymbol()}"
    
    return True, "Molecule is a nitrohydrocarbon: It is a hydrocarbon with at least one nitro group replacing an H."

    
# The following block is for example testing. In actual use this function could be imported into a larger program.
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
        # Expected false positives:
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