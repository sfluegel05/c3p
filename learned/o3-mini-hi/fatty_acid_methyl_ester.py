"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: Fatty acid methyl ester
A fatty acid methyl ester is defined as a fatty acid ester that is the carboxylic ester obtained
by the formal condensation of a fatty acid with methanol.
The algorithm first finds the methyl ester substructure [CX3](=O)O[CH3] and then examines 
the acyl portion (the R in R-C(=O)OCH3). For an R-group to be considered a fatty acid chain, 
we require that it consists primarily of carbon atoms and that if it is long (≥5 carbons) it is
sufficiently flexible (i.e. contains many rotatable bonds). For short chains (<5 carbons) we use 
a minimum carbon count criterion.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    It first looks for the methyl ester group [CX3](=O)O[CH3]. Then it identifies the acyl chain
    attached to the carbonyl carbon (excluding the carbonyl O) and computes its characteristics:
      - number of carbon atoms in the chain (chain_length)
      - number of rotatable bonds in that chain (acyl_rotatable).
      
    For chains with at least 5 carbon atoms, we require a minimum flexibility 
    (acyl_rotatable >= chain_length*0.5, rounded to an integer) to qualify as a fatty acid. 
    For chains shorter than 5 carbons, we accept a minimum chain length of 3 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the classification criteria are met, False otherwise
        str: Reason for classification
    """
    # Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS for a methyl ester group: C(=O)OCH3.
    ester_smarts = "[CX3](=O)O[CH3]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error defining ester SMARTS"

    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No methyl ester group found in the molecule"

    # Helper function: Given a starting carbon atom (the acyl chain) and an atom to exclude,
    # perform a breadth-first search to collect indices of connected carbon atoms.
    def get_acyl_chain_atoms(start_atom, exclude_idx):
        chain_atoms = set()
        queue = [start_atom.GetIdx()]
        while queue:
            idx = queue.pop(0)
            if idx in chain_atoms:
                continue
            chain_atoms.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            # Traverse only if the atom is carbon (atomic num == 6).
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == exclude_idx:
                    continue
                # Only follow carbon neighbors (allowing unsaturation, etc.)
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in chain_atoms:
                    queue.append(nbr.GetIdx())
        return chain_atoms

    # Helper function: Count rotatable bonds in the acyl chain.
    # We consider a bond rotatable if it is a single bond (bond order == 1), not in a ring,
    # and connects two heavy atoms that are not terminal (i.e. having more than one neighbor in the chain).
    def count_rotatable_bonds_in_chain(chain_atom_indices):
        rot_bonds = 0
        # Create a set for quick lookup.
        chain_set = set(chain_atom_indices)
        for bond in mol.GetBonds():
            # Check if both atoms belong to the acyl chain.
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetIdx() in chain_set and a2.GetIdx() in chain_set:
                # Check for single bond and not in ring.
                if bond.GetBondTypeAsDouble() == 1.0 and not bond.IsInRing():
                    # Check that neither end is terminal within the chain. Count neighbors that are in chain.
                    a1_neighbors_in_chain = sum(1 for nbr in a1.GetNeighbors() if nbr.GetIdx() in chain_set)
                    a2_neighbors_in_chain = sum(1 for nbr in a2.GetNeighbors() if nbr.GetIdx() in chain_set)
                    if a1_neighbors_in_chain > 1 and a2_neighbors_in_chain > 1:
                        rot_bonds += 1
        return rot_bonds

    # Process each found methyl ester match.
    # Each match is a tuple of atom indices corresponding to the SMARTS:
    # (carbonyl carbon, carbonyl oxygen, ester oxygen, methyl carbon)
    for match in ester_matches:
        carbonyl_idx, carbonylO_idx, esterO_idx, methyl_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the acyl neighbor: find the neighbor of the carbonyl that is not the carbonyl O.
        acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == carbonylO_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_neighbor = nbr
                break
        if acyl_neighbor is None:
            continue  # No acyl chain here; try next match.
        
        # Collect all connected carbon atoms that form the acyl chain.
        acyl_atoms = get_acyl_chain_atoms(acyl_neighbor, carbonyl_idx)
        chain_length = len(acyl_atoms)
        
        # For a fatty acid, even very short ones must have at least 3 carbons.
        if chain_length < 3:
            continue
        
        # Count rotatable bonds in the acyl chain.
        acyl_rotatable = count_rotatable_bonds_in_chain(acyl_atoms)
        
        # Apply flexibility criteria:
        # For chains with 5 or more carbons, we require a reasonable number of rotatable bonds.
        # (A loose heuristic is: rotatable bonds should be at least half of chain_length, rounded down.)
        if chain_length >= 5:
            required_rot_bonds = chain_length // 2
            if acyl_rotatable < required_rot_bonds:
                continue
        
        # We can also check some overall molecule descriptors for reasonability.
        n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        reason = (f"Found a methyl ester group with an acyl chain of {chain_length} carbon(s); "
                  f"chain rotatable bonds: {acyl_rotatable}, overall rotatable bonds: {n_rotatable}, "
                  f"molecular weight: {mol_wt:.1f} Da")
        return True, reason

    return False, "Molecule contains a methyl ester group but the acyl (fatty acid) part does not meet the criteria"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "O=C(OC)CCCCCCCCCCCC",  # methyl tridecanoate, expected True
        "COC(=O)[C@@H]1O[C@@]1(C)CC\\C=C(/C)CC[C@H]1OC1(C)C",  # juvenile hormone III skipped bisepoxide, expected True
        "COC(=O)CCC(=O)CN",  # methyl 5-aminolevulinate, expected False
        "BrC1=C(O)C2=C([C@@H]3[C@](CC[C@@H]3C(C)C)(C)C([C@]2(O)C(=O)OC)=O)C=C1C"  # false positive gibberellin A3 methyl ester, expected False
    ]
    for sm in test_smiles:
        valid, msg = is_fatty_acid_methyl_ester(sm)
        print(f"SMILES: {sm}\nResult: {valid}\nReason: {msg}\n")