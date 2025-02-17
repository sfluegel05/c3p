"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17α-hydroxy steroid
Definition: The α-stereoisomer of 17-hydroxy steroid.
This revised code employs a more stringent heuristic:
  1. It parses the input SMILES.
  2. It computes the molecule’s Murcko scaffold (the “core”).
  3. It checks that the scaffold contains exactly 4 fused rings – three six‐membered rings and one five‐membered ring.
  4. It further requires that the scaffold contains exactly 17 carbon atoms (the classical steroid nucleus).
  5. Finally, it looks for at least one hydroxyl group (–OH) that is attached to a chiral carbon residing in the five‐membered ring (D–ring).
If all tests are passed the molecule is heuristically classified as a 17α-hydroxy steroid.
Note: Due to practical limitations a full assignment of the C17 position and its α–configuration is not performed.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17α-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a 17α-hydroxy steroid, False otherwise
        str: Explanation of the classification decision
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain the Murcko scaffold (the core structure) of the molecule.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not obtain Murcko scaffold"
    
    # Get the ring information of the scaffold.
    rings = scaffold.GetRingInfo().AtomRings()
    if len(rings) != 4:
        return False, ("Molecule's core does not have exactly 4 fused rings typical of the steroid nucleus. "
                       f"Found {len(rings)} rings in the scaffold.")
    
    # Check that the scaffold has exactly three six‐membered rings and one five‐membered ring.
    ring_sizes = [len(r) for r in rings]
    if ring_sizes.count(6) != 3 or ring_sizes.count(5) != 1:
        return False, ("Molecule's core rings do not match expected sizes (3 six‐membered and 1 five‐membered). "
                       f"Found ring sizes: {ring_sizes}")
    
    # Check that the scaffold contains exactly 17 carbon atoms (a common feature of the steroid nucleus).
    carbon_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 17:
        return False, (f"Molecule's scaffold does not have 17 carbons typical of steroids (found {carbon_count}).")
    
    # To ensure that the hydroxyl is at the appropriate position we now require that an -OH is
    # attached to a chiral carbon that belongs to the five-membered ring in the scaffold.
    # First, get the indices (in the scaffold) of the atoms in the five-membered ring.
    five_ring_indices_scaffold = set()
    for r in rings:
        if len(r) == 5:
            five_ring_indices_scaffold.update(r)
    
    # Map the scaffold atoms back to the original molecule.
    scaffold_matches = mol.GetSubstructMatches(scaffold)
    if not scaffold_matches:
        return False, "Unable to map the calculated scaffold to the molecule"
    # Use the first mapping (a tuple where the i-th scaffold atom corresponds to mol atom index mapping[i])
    mapping = scaffold_matches[0]
    # Build a set of molecule atom indices that belong to the five-membered ring of the scaffold.
    five_ring_mol_indices = set(mapping[i] for i in five_ring_indices_scaffold)
    
    # Now search for hydroxyl groups (-OH) in the molecule.
    oh_smarts = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_smarts)
    if not oh_matches:
        return False, "No hydroxyl (–OH) group found"
    
    candidate_found = False
    # Check each hydroxyl group: the oxygen must be connected to a carbon that (a) is part of the five‐membered ring and (b) is chiral.
    for match in oh_matches:
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        for nb in o_atom.GetNeighbors():
            if nb.GetAtomicNum() == 6 and nb.GetIdx() in five_ring_mol_indices:
                # Check if the carbon has an assigned chirality.
                if nb.HasProp('_CIPCode'):
                    candidate_found = True
                    break
        if candidate_found:
            break
    if not candidate_found:
        return False, ("No hydroxyl group was found at a chiral carbon in the expected five-membered (D-ring) "
                       "of the steroid nucleus.")
    
    # All heuristics passed; assume the molecule is a 17α-hydroxy steroid.
    return True, ("Molecule’s Murcko scaffold shows a classical steroid nucleus (3 six‐membered and 1 five‐membered rings, "
                  "17 carbons) and it has a hydroxyl group attached to a chiral carbon in the D‐ring; "
                  "heuristically classified as a 17α-hydroxy steroid.")

# Example usage (uncomment the following lines to test):
# test_smiles = "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO"  # methyl prednisolone-16alpha-carboxylate
# result, reason = is_17alpha_hydroxy_steroid(test_smiles)
# print(result, reason)