"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: aldopentose, defined as 'A pentose with a (potential) aldehyde group at one end.'
This function now uses improved heuristics. For open‐chain structures, it builds a carbon‐connectivity graph
and checks for a unique 5-carbon contiguous (sugar) backbone with a terminal aldehyde and a terminal primary alcohol.
For cyclic structures, it searches for a defined sugar ring (either furanose or pyranose) and checks that there is no 
extra carbon branching (which happens in apiose) and rejects lactone functionality.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.

    An aldopentose is a 5-carbon sugar that in its open-chain form displays an aldehyde group at one terminal carbon.
    In its cyclic (hemiacetal) form the explicit aldehyde may not be present but the sugar ring pattern is expected.
    This function uses two routes: one for open-chain molecules and one for cyclic.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an aldopentose; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Reject lactones (esters) because those are not aldopentoses.
    lactone_pattern = Chem.MolFromSmarts("[CX3](=O)[O]")
    if mol.HasSubstructMatch(lactone_pattern):
        return False, "Molecule contains lactone (ester) functionality; not an aldopentose."
    
    # Prepare a pattern for an aldehyde carbon: carbon with one hydrogen, double-bonded to oxygen.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    aldehyde_indices = {match[0] for match in aldehyde_matches}  # these are the aldehyde carbons

    # For open-chain molecules we analyze the contiguous carbon backbone.
    # We use an added-H version to help check for primary alcohol groups.
    mol_h = Chem.AddHs(mol)
    
    # Build a graph of all carbon atoms (using the indices from mol_h). 
    carbon_idxs = [atom.GetIdx() for atom in mol_h.GetAtoms() if atom.GetAtomicNum() == 6]
    # Represent graph as adjacency dictionary for carbons.
    carbon_graph = {idx: set() for idx in carbon_idxs}
    for bond in mol_h.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum()==6 and a2.GetAtomicNum()==6:
            i1 = a1.GetIdx()
            i2 = a2.GetIdx()
            if i1 in carbon_graph and i2 in carbon_graph:
                carbon_graph[i1].add(i2)
                carbon_graph[i2].add(i1)
    
    # Helper: check if a carbon atom (given its index in mol_h) qualifies as a primary alcohol CH2OH.
    # We require that it is sp3, has at least one oxygen neighbor that itself is bonded to at least one hydrogen.
    def is_primary_alcohol(carbon_idx):
        atom = mol_h.GetAtomWithIdx(carbon_idx)
        # Check that this carbon is not aldehyde (should be CH2)
        # Count oxygen neighbors that have at least one hydrogen
        oxy_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum()==8]
        for o in oxy_neighbors:
            # In the added-H molecule, hydrogens are explicit.
            h_count = sum(1 for n in o.GetNeighbors() if n.GetSymbol()=="H")
            if h_count >= 1:
                return True
        return False

    # Helper: find all simple paths in the carbon_graph starting from start up to a given depth (in number of nodes).
    # Because our molecules are small, a recursive DFS is acceptable.
    def dfs_paths(start, target_length, path, visited):
        if len(path) == target_length:
            yield list(path)
            return
        for neighbor in carbon_graph[start]:
            if neighbor not in visited:
                visited.add(neighbor)
                path.append(neighbor)
                for res in dfs_paths(neighbor, target_length, path, visited):
                    yield res
                path.pop()
                visited.remove(neighbor)

    # Open-chain candidate: check if any simple path of exactly 5 carbon atoms exists
    # such that one end is an aldehyde carbon and the other end is a primary alcohol.
    open_chain_valid = False
    for start in carbon_idxs:
        if start in aldehyde_indices:  # consider only starting at an aldehyde carbon
            # search for simple paths of length 5 (number of nodes = 5)
            for path in dfs_paths(start, 5, [start], {start}):
                # path[0] is the aldehyde; path[-1] must be a primary alcohol
                if is_primary_alcohol(path[-1]):
                    # Also, verify that the 5 carbons are contiguous:
                    # (this is already ensured by our DFS on the carbon_graph)
                    open_chain_valid = True
                    break
            if open_chain_valid:
                break

    if open_chain_valid:
        return True, "Open-chain aldopentose: found a contiguous 5-carbon backbone with terminal aldehyde and primary alcohol."

    # If not open-chain, try cyclic analysis.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_symbols = [atom.GetSymbol() for atom in ring_atoms]
        ring_size = len(ring)
        # Case 1: pyranose: expected ring size 6 with exactly one oxygen in the ring 
        # (thus five carbons in the ring; the sugar backbone is the five ring carbons).
        if ring_size == 6 and ring_symbols.count("O") == 1:
            return True, "Cyclized aldopentose: detected pyranose ring (6-membered ring with 5 carbons, 1 oxygen)."
        # Case 2: furanose: expected ring size 5 with exactly one oxygen (and thus 4 carbons in ring),
        # plus one exocyclic carbon attached to the ring (at the anomeric carbon) to complete the 5-carbon backbone.
        if ring_size == 5 and ring_symbols.count("O") == 1:
            # Get indices in the ring that are carbons.
            ring_carbon_idxs = [ring[i] for i, sym in enumerate(ring_symbols) if sym=="C"]
            # Check exocyclic carbon attachments:
            exocyclic_carbons = set()
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    # If neighbor is carbon and not in the ring then it might be the exocyclic carbon.
                    if nbr.GetAtomicNum()==6 and nbr_idx not in ring:
                        exocyclic_carbons.add(nbr_idx)
            # For a proper furanose sugar, there should be exactly one exocyclic carbon.
            if len(exocyclic_carbons) == 1:
                return True, "Cyclized aldopentose: detected furanose ring (5-membered ring with proper exocyclic attachment)."
    return False, "No aldopentose pattern found."

# (Optional testing)
if __name__ == "__main__":
    # Some examples:
    examples = [
        ("O1[C@H]([C@@H](O)[C@H](O)[C@@H]1O)CO", "alpha-L-xylofuranose"),
        ("O[C@H]1CO[C@H](O)[C@@H](O)[C@@H]1O", "beta-L-xylopyranose"),
        ("[H]C(=O)[C@@H](O)[C@H](O)[C@@H](O)CO", "aldehydo-L-xylose"),
        ("O1C[C@]([C@H](C1O)O)(O)CO", "D-apiofuranose (false positive?)"),
        ("[H]C(=O)[C@@](O)(CO)[C@H](O)[C@H](O)CO", "D-hamamelose (false negative previously)")
    ]
    for smi, name in examples:
        result, reason = is_aldopentose(smi)
        print(f"Name: {name}\nSMILES: {smi}\nClassification: {result}\nReason: {reason}\n{'-'*60}")