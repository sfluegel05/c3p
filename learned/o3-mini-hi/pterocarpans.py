"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: Pterocarpans – Members of the class of benzofurochromene with a
6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton (the 3,4-dihydro derivatives of coumestans).
We now use several heuristics:
  1. Extract the Bemis–Murcko scaffold.
  2. Require that the scaffold accounts for a large fraction of the molecule’s heavy atoms.
  3. Accept scaffolds with a modest number of fused rings (3–5) so that additional fused sub‐rings are allowed.
  4. Require the scaffold contains at least one oxygen.
  5. Require a benzofuran substructure (SMARTS: "c1ccc2occc2c1") as pterocarpans belong to the benzofurochromene family.
  6. Require that at least one ring is “dihydro” (contains a non‐aromatic carbon) to mimic the 3,4-dihydro character.
This heuristic is not perfect, but we hope it better rejects false positives and identifies known pterocarpans.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is likely a pterocarpan based on its SMILES string.
    
    Approach:
      1. Parse the molecule.
      2. Extract its Bemis–Murcko scaffold.
      3. Check that the scaffold makes up a large fraction of the molecule.
      4. Accept only scaffolds with between 3 and 5 rings (to allow extra fused rings from derivatization).
      5. Verify the scaffold contains at least one oxygen atom.
      6. Verify the scaffold contains a benzofuran substructure.
      7. Ensure that at least one ring has one or more saturated (non‐aromatic) carbons to capture the 3,4-dihydro feature.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a pterocarpan, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extract the Bemis-Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error extracting scaffold: {e}"
    if scaffold is None:
        return False, "Could not extract scaffold"

    # Criterion 1: The scaffold should represent a large fraction of the molecule.
    mol_heavy = mol.GetNumHeavyAtoms()
    scaffold_heavy = scaffold.GetNumHeavyAtoms()
    ratio = scaffold_heavy / mol_heavy if mol_heavy > 0 else 0
    if ratio < 0.6:
        return False, f"Scaffold heavy atom ratio too low ({ratio:.2f}); likely not a pterocarpan core"

    # Criterion 2: Check that the scaffold has a number of rings consistent with a pterocarpan core.
    ri = scaffold.GetRingInfo()
    atom_rings = ri.AtomRings()  # Each ring as a tuple of atom indices.
    n_rings = len(atom_rings)
    if n_rings < 3 or n_rings > 5:
        return False, f"Scaffold has {n_rings} ring(s); expected between 3 and 5 fused rings for a pterocarpan core"

    # Criterion 3: Check that the scaffold contains at least one oxygen atom.
    oxygens = [atom for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 8]
    if not oxygens:
        return False, "Scaffold does not contain any oxygen atoms"

    # Criterion 4: Check for a benzofuran substructure.
    # The SMARTS "c1ccc2occc2c1" represents a benzofuran system.
    benzofuran_smarts = "c1ccc2occc2c1"
    benzofuran = Chem.MolFromSmarts(benzofuran_smarts)
    if not scaffold.HasSubstructMatch(benzofuran):
        return False, "Scaffold does not contain a benzofuran substructure required for pterocarpans"

    # Criterion 5: Check for a dihydro (reduced) ring: at least one ring should have a non‐aromatic (saturated) carbon.
    ring_has_saturated = False
    for ring in atom_rings:
        for idx in ring:
            atom = scaffold.GetAtomWithIdx(idx)
            # If a carbon atom is not aromatic, mark the ring as reduced.
            if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
                ring_has_saturated = True
                break
        if ring_has_saturated:
            break
    if not ring_has_saturated:
        return False, "No dihydro (reduced) ring found in the scaffold"

    # Passed all criteria.
    return True, ("Molecule meets heuristic criteria for pterocarpans: scaffold covers a high fraction of the molecule, "
                  "has 3–5 fused rings containing oxygen, shows a benzofuran substructure, and has a dihydro ring.")

    
# Example usage (for testing):
if __name__ == "__main__":
    test_examples = [
        # True positives (pterocarpans)
        ("O1C2C(C3=C1C=C(O)C=C3)COC4=C2C=CC(O)=C4O", "3,4,9-Trihydroxypterocarpan"),
        ("O1C2C(C=3C1=C(OC)C(OC)=C(O)C3)COC=4C2=CC(O)=C(OC)C4", "2,8-Dihydroxy-3,9,10-trimethoxypterocarpan"),
        ("COc1c(CC=C(C)C)c(O)cc2OC[C@@H]3[C@@H](Oc4cc(O)ccc34)c12", "edudiol"),
        ("O1C2C(C3=C1C(OC)=C(OC)C=C3)COC4=C2C=CC(OC)=C4O", "Melilotocarpan C"),
        # Some additional examples that were problematic:
        ("[H][C@@]12COc3cc(O)c(CC=C(C)C)cc3[C@]1([H])Oc1c(CC=C(C)C)c(OC)c(O)cc21", "lespeflorin G2"),
        ("O1C2C(C=3C1=CC(OC)=CC3O)COC4=C2C=CC(O)=C4", "Nissicarpin"),
        # False positive example (non-pterocarpan)
        ("COC1=CC=C(C=C1)C(=O)CSC2=NN=C(N2CC=C)C3=CC4=C(O3)C(=CC=C4)OC", "Non-pterocarpan example")
    ]
    
    for smi, name in test_examples:
        res, reason = is_pterocarpans(smi)
        print(f"Name: {name}")
        print(f"SMILES: {smi}")
        print(f"Result: {res}")
        print(f"Reason: {reason}\n")