"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: Pterocarpans – Members of the class of benzofurochromene with a
6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton (the 3,4-dihydro derivatives 
of coumestans) and its substituted derivatives.
Heuristics used here:
  1. Extract the Bemis–Murcko scaffold.
  2. Require that the scaffold represents a significant fraction of the molecule (heavy atom ratio ≥ 0.5).
  3. Accept only scaffolds with 3–5 rings (to allow extra fused sub‐rings due to derivatization).
  4. Require that the scaffold contains at least one oxygen atom.
  5. Require that either the scaffold or the entire molecule contains a benzofuran substructure 
     (SMARTS: "c1ccc2occc2c1").
  6. Require that at least one ring in the scaffold has a saturated (non‐aromatic) carbon to capture 
     the 3,4-dihydro (reduced) character.
Note: This heuristic is not perfect but attempts to balance false positives and negatives.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is likely a pterocarpan based on its SMILES string using a set of heuristics.
    
    Heuristics:
      1. Extract the Bemis–Murcko scaffold and compute the heavy atom ratio (scaffold/molecule) which must be ≥ 0.5.
      2. The scaffold must have 3–5 rings.
      3. The scaffold must contain at least one oxygen atom.
      4. Either the scaffold or the full molecule must match the benzofuran substructure (SMARTS: "c1ccc2occc2c1").
      5. At least one ring in the scaffold must contain a saturated (non‐aromatic) carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a pterocarpan, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extract Bemis–Murcko scaffold
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error extracting scaffold: {e}"
    if scaffold is None:
        return False, "Could not extract scaffold"
    
    # Criterion 1: Check heavy atom ratio (relaxed to 0.5)
    mol_heavy = mol.GetNumHeavyAtoms()
    scaffold_heavy = scaffold.GetNumHeavyAtoms()
    if mol_heavy == 0:
        return False, "Molecule has no heavy atoms"
    ratio = scaffold_heavy / mol_heavy
    if ratio < 0.5:
        return False, f"Scaffold heavy atom ratio too low ({ratio:.2f}); likely not a pterocarpan core"
    
    # Criterion 2: Check for 3-5 rings in the scaffold.
    ri = scaffold.GetRingInfo()
    atom_rings = ri.AtomRings()  # Each ring is a tuple of atom indices.
    n_rings = len(atom_rings)
    if n_rings < 3 or n_rings > 5:
        return False, f"Scaffold has {n_rings} ring(s); expected between 3 and 5 fused rings for a pterocarpan core"
    
    # Criterion 3: Scaffold must have at least one oxygen.
    oxygens = [atom for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 8]
    if not oxygens:
        return False, "Scaffold does not contain any oxygen atoms"
    
    # Criterion 4: Check for benzofuran substructure.
    benzofuran_smarts = "c1ccc2occc2c1"
    benzofuran = Chem.MolFromSmarts(benzofuran_smarts)
    # We require that either the scaffold OR the original molecule has this substructure.
    if not (scaffold.HasSubstructMatch(benzofuran) or mol.HasSubstructMatch(benzofuran)):
        return False, "Neither scaffold nor full molecule contains a benzofuran substructure required for pterocarpans"
    
    # Criterion 5: Require at least one ring in the scaffold has a saturated (non‐aromatic) carbon.
    ring_has_saturated = False
    for ring in atom_rings:
        # For each ring, check if there is at least one carbon that is not aromatic.
        for idx in ring:
            atom = scaffold.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
                ring_has_saturated = True
                break
        if ring_has_saturated:
            break
    if not ring_has_saturated:
        return False, "No dihydro (reduced) ring found in the scaffold"
    
    # If all criteria are passed.
    return True, ("Molecule meets heuristic criteria for pterocarpans: scaffold covers a high fraction of the molecule, "
                  "has 3–5 fused rings with an oxygen, contains a benzofuran substructure (in scaffold or molecule), "
                  "and has a dihydro ring (with at least one saturated carbon).")

    
# For testing purposes:
if __name__ == "__main__":
    test_examples = [
        # True positives (pterocarpans)
        ("O1C2C(O)(C3=C1C=C(O)C=C3)COC4=C2C=CC(O)=C4O", "3,4,9-Trihydroxypterocarpan"),
        ("O1C2C(C=3C1=C(OC)C(OC)=C(O)C3)COC=4C2=CC(O)=C(OC)C4", "2,8-Dihydroxy-3,9,10-trimethoxypterocarpan"),
        ("COc1c(CC=C(C)C)c(O)cc2OC[C@@H]3[C@@H](Oc4cc(O)ccc34)c12", "edudiol"),
        ("O1C2C(C3=C1C(OC)=C(OC)C=C3)COC4=C2C=CC(OC)=C4O", "Melilotocarpan C"),
        ("[H][C@@]12COc3cc(O)c(CC=C(C)C)cc3[C@]1([H])Oc1c(CC=C(C)C)c(OC)c(O)cc21", "lespeflorin G2"),
        ("O1C2C(C=3C1=CC(OC)=CC3O)COC4=C2C=CC(O)=C4", "Nissicarpin"),
        # A false positive example (non-pterocarpan but may share substructures)
        ("COC1=CC=C(C=C1)C(=O)CSC2=NN=C(N2CC=C)C3=CC4=C(O3)C(=CC=C4)OC", "Non-pterocarpan example")
    ]
    
    for smi, name in test_examples:
        res, reason = is_pterocarpans(smi)
        print(f"Name: {name}")
        print(f"SMILES: {smi}")
        print(f"Result: {res}")
        print(f"Reason: {reason}\n")