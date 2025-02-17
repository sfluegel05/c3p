"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: Pterocarpans – Members of the class of benzofurochromene with a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.
They are 3,4-dihydro derivatives of coumestans and generally resemble isoflavanoids produced by plants.
We attempt a classification by extracting the Murcko scaffold and checking for:
  (1) Exactly three fused rings (the pterocarpan core is a tricyclic system),
  (2) The scaffold contains at least one oxygen atom,
  (3) At least one of the rings has one or more non-aromatic (saturated) carbon atoms (reflecting the dihydro nature).
This is a heuristic and may mis‐classify unusual derivatives.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    
    Our approach is:
      1. Parse the SMILES into a molecule.
      2. Extract its Bemis-Murcko scaffold.
      3. Count the number of rings in the scaffold.
      4. Check that the scaffold contains at least one oxygen.
      5. Identify that (at least) one of the rings is not fully aromatic (i.e. has a saturated carbon)
         to capture the 3,4-dihydro (reduced) character.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the scaffold meets our pterocarpan criteria, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the Murcko scaffold (the “core” ring system after removal of sidechains)
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error extracting scaffold: {e}"
    
    if scaffold is None:
        return False, "Could not extract scaffold"
    
    # Compute ring info from the scaffold
    ri = scaffold.GetRingInfo()
    atom_rings = ri.AtomRings()  # tuple of atom-index tuples for each ring
    num_rings = len(atom_rings)
    
    # We expect the pterocarpan core to have three fused rings.
    if num_rings != 3:
        return False, f"Scaffold has {num_rings} ring(s); expected 3 fused rings for pterocarpan core"
    
    # Check that the scaffold contains at least one oxygen atom.
    oxygens = [atom for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 8]
    if not oxygens:
        return False, "Scaffold does not contain any oxygen atoms"
    
    # For each ring in the scaffold, check if it is fully aromatic.
    # We expect a 3,4-dihydro derivative: one ring should have at least one carbon that is not aromatic.
    ring_has_saturated = [False] * num_rings
    for i, ring in enumerate(atom_rings):
        for idx in ring:
            atom = scaffold.GetAtomWithIdx(idx)
            # Consider carbons that are not aromatic as “saturated”
            if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()):
                ring_has_saturated[i] = True
                break

    if not any(ring_has_saturated):
        return False, "All rings in the scaffold are aromatic; no dihydro (saturated) ring found"
    
    # If we have three fused rings, at least one oxygen, and at least one reduced ring, then we classify as pterocarpan.
    return True, "Scaffold has 3 fused rings (benzofurochromene core) with at least one dihydro (reduced) ring and oxygen."

    
# Example usage (for testing):
if __name__ == "__main__":
    test_examples = [
        # True positives (expected to be pterocarpans)
        ("O1C2C(C3=C1C=C(O)C=C3)COC4=C2C=CC(O)=C4O", "3,4,9-Trihydroxypterocarpan"),
        ("O1C2C(C=3C1=C(OC)C(OC)=C(O)C3)COC=4C2=CC(O)=C(OC)C4", "2,8-Dihydroxy-3,9,10-trimethoxypterocarpan"),
        ("COc1c(CC=C(C)C)c(O)cc2OC[C@@H]3[C@@H](Oc4cc(O)ccc34)c12", "edudiol"),
        ("O1C2C(C3=C1C(OC)=C(OC)C=C3)COC4=C2C=CC(OC)=C4O", "Melilotocarpan C"),
        # A false example (likely not a pterocarpan)
        ("COC1=CC=C(C=C1)C(=O)CSC2=NN=C(N2CC=C)C3=CC4=C(O3)C(=CC=C4)OC", "Non-pterocarpan example")
    ]
    
    for smi, name in test_examples:
        res, reason = is_pterocarpans(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n")