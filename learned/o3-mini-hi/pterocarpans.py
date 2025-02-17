"""
Classifies: CHEBI:26377 pterocarpans
"""
#!/usr/bin/env python
"""
Classifies: Pterocarpans
Definition: Members of the class of benzofurochromene with a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton 
           (i.e. the 3,4-dihydro derivatives of coumestans) and its substituted derivatives.
Heuristic improvements in this version:
  1. Extract the Bemis–Murcko scaffold and require that it represents at least 50% of the molecule’s heavy atoms.
  2. Require that the scaffold has at least 3 rings.
  3. Require that the scaffold contains at least one oxygen atom.
  4. Require that at least 2 of the rings in the scaffold are fully aromatic (fused aromatic core).
  5. Require that at least one ring in the scaffold is “dihydro” (contains at least one carbon that is non‐aromatic) 
     to capture the reduced (3,4-dihydro) character.
This heuristic tries to balance capturing the structural essence of pterocarpans while avoiding molecules 
that merely share one or two substructures.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is likely a pterocarpan based on its SMILES string using improved heuristics.
    
    Heuristics:
      1. Extract the Bemis–Murcko scaffold and compute the heavy atom ratio (scaffold/molecule) which must be ≥ 0.50.
      2. The scaffold must have at least 3 rings.
      3. The scaffold must contain at least one oxygen atom.
      4. The scaffold must contain at least 2 fully aromatic rings, representing the fused aromatic portion.
      5. At least one ring in the scaffold must have at least one saturated (non‐aromatic) carbon to reflect its
         dihydro (reduced) character.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a pterocarpan, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extract Bemis–Murcko scaffold from molecule
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error extracting scaffold: {e}"
    if scaffold is None:
        return False, "Could not extract scaffold from molecule"
    
    # Criterion 1: Check the heavy atom ratio (scaffold vs full molecule)
    mol_heavy = mol.GetNumHeavyAtoms()
    scaffold_heavy = scaffold.GetNumHeavyAtoms()
    if mol_heavy == 0:
        return False, "Molecule has no heavy atoms"
    ratio = scaffold_heavy / mol_heavy
    if ratio < 0.50:
        return False, f"Scaffold heavy atom ratio too low ({ratio:.2f}); core may not be representative"
    
    # Criterion 2: The scaffold must have at least 3 rings.
    ring_info = scaffold.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # list of tuples (each ring: indices of atoms)
    n_rings = len(atom_rings)
    if n_rings < 3:
        return False, f"Scaffold has {n_rings} ring(s); expected at least 3 fused rings for a pterocarpan core"
    
    # Criterion 3: Scaffold must contain at least one oxygen atom.
    if not any(atom.GetAtomicNum() == 8 for atom in scaffold.GetAtoms()):
        return False, "Scaffold does not contain any oxygen atoms"
    
    # Criterion 4: Require that the scaffold has at least two fully aromatic rings.
    aromatic_ring_count = 0
    for ring in atom_rings:
        # Check if every atom in the ring is aromatic.
        if all(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    if aromatic_ring_count < 2:
        return False, f"Scaffold has only {aromatic_ring_count} fully aromatic ring(s); expected at least 2 for the fused aromatic core"
    
    # Criterion 5: Check for presence of a "dihydro" (reduced) ring.
    # That is, at least one ring in the scaffold must contain at least one saturated (non‐aromatic) carbon.
    has_dihydro = False
    for ring in atom_rings:
        for idx in ring:
            atom = scaffold.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()):
                # This carbon is not part of an aromatic system, indicative of a partially saturated ring.
                has_dihydro = True
                break
        if has_dihydro:
            break
    if not has_dihydro:
        return False, "No dihydro (reduced) ring found in the scaffold"
    
    # All criteria passed—molecule meets our heuristic criteria for pterocarpans.
    return True, ("Molecule meets heuristic criteria for pterocarpans: scaffold covers a sufficient fraction of the molecule, "
                  "has at least 3 rings (with at least 2 fully aromatic ones), contains an oxygen, and shows dihydro (reduced) features.")

# Optional testing block
if __name__ == "__main__":
    test_examples = [
        # True positive examples:
        ("O1C2C(O)(C3=C1C=C(O)C=C3)COC4=C2C=CC(O)=C4O", "3,4,9-Trihydroxypterocarpan"),
        ("O1C2C(C=3C1=C(OC)C(OC)=C(O)C3)COC=4C2=CC(O)=C(OC)C4", "2,8-Dihydroxy-3,9,10-trimethoxypterocarpan"),
        ("COc1c(CC=C(C)C)c(O)cc2OC[C@@H]3[C@@H](Oc4cc(O)ccc34)c12", "edudiol"),
        ("O1C2C(C3=C1C(OC)=C(OC)C=C3)COC4=C2C=CC(OC)=C4O", "Melilotocarpan C"),
        # False positive example (a molecule that may share parts but is not a pterocarpan):
        ("O1C2=C3C(=C(C)C=C2C(=C1)CO)CC[C@H]([C@@H]3C)O", "Citreobenzofuran B")
    ]
    
    for smi, name in test_examples:
        res, reason = is_pterocarpans(smi)
        print(f"Name: {name}")
        print(f"SMILES: {smi}")
        print(f"Result: {res}")
        print(f"Reason: {reason}\n")