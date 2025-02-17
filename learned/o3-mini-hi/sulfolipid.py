"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: A sulfolipid, defined as:
'A compound containing a sulfonic acid residue joined by a carbon–sulfur bond to a lipid.'
In practice, many sulfolipids (e.g. sulfated glycolipids or sphingolipids) are large molecules
with a polar headgroup (frequently a sugar or sphingosine ring system) plus long aliphatic chains.
This implementation uses several heuristics:
 – Molecular weight must be above 300 Da.
 – The molecule must contain at least one ring.
 – A sulfur atom is searched that exhibits at least two S=O bonds and one S–O (or S–OH) bond.
 – The sulfur must be connected (either directly or via an oxygen) to a carbon that is part of
   a long (≥12 carbon) aliphatic chain.
If all these conditions are met, the molecule is classified as a sulfolipid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is defined as a compound containing a sulfonic acid residue (S(=O)(=O)[O] or S(=O)(=O)O)
    attached (directly or indirectly via a single oxygen) to a long aliphatic (lipid) chain.
    Additionally, most known sulfolipids possess extra polar headgroup features such as rings.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a sulfolipid, False otherwise.
        str: Reason for the classification result.
    """
    # First, parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic: require a minimum molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a sulfolipid (mol wt = {mol_wt:.1f})"
    
    # Heuristic: many sulfolipids have a polar headgroup (often a sugar or sphingoid)
    # so require at least one ring.
    if len(mol.GetRingInfo().AtomRings()) == 0:
        return False, "No ring structures found (expected polar headgroup feature)"
    
    # Define a DFS routine to compute the length of a contiguous aliphatic chain.
    def dfs_chain(atom, visited):
        count = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in visited:
                continue
            # Only count aliphatic carbons: atomic number 6, sp3, non-aromatic, not in a ring.
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() and not nbr.IsInRing():
                new_visited = visited.copy()
                new_visited.add(nbr.GetIdx())
                candidate = 1 + dfs_chain(nbr, new_visited)
                if candidate > count:
                    count = candidate
        return count

    # Now search for candidate sulfur atoms.
    # We look for a sulfur that has at least two double bonds to oxygen and one single bond to oxygen.
    # (Sometimes the single-bond oxygen can be deprotonated so the charge may vary.)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # Not sulfur.
        neighbors = atom.GetNeighbors()
        
        o_double_count = 0
        o_single_count = 0
        attached_atoms = []  # list of neighbor atoms that are not oxygen (or, for bridging O, we will check later)
        
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # Use GetBondTypeAsDouble() to distinguish double vs single bonds.
                if abs(bond.GetBondTypeAsDouble() - 2.0) < 0.01:
                    o_double_count += 1
                elif abs(bond.GetBondTypeAsDouble() - 1.0) < 0.01:
                    o_single_count += 1
            else:
                attached_atoms.append(nbr)
        # We require at least two double bonds and one single bond to oxygen.
        if o_double_count < 2 or o_single_count < 1:
            continue  # This S is not clearly in a sulfonic acid environment.
        
        # Now, among the attached atoms (i.e. non-oxygen neighbors), check for a carbon atom.
        # Also, sometimes the S may be attached indirectly via an oxygen (C-O-S pattern).
        candidate_c_atoms = []
        for nbr in attached_atoms:
            if nbr.GetAtomicNum() == 6:
                candidate_c_atoms.append(nbr)
        
        # Also, explore if S is connected to an oxygen that in turn is connected to a carbon.
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 8:
                # (Only follow if the oxygen is single-bonded to S.)
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None or abs(bond.GetBondTypeAsDouble()-1.0) > 0.01:
                    continue
                for nn in nbr.GetNeighbors():
                    if nn.GetAtomicNum() == 6 and nn.GetIdx() != atom.GetIdx():
                        candidate_c_atoms.append(nn)
        
        # If no candidate carbon is found, continue to the next sulfur.
        if not candidate_c_atoms:
            continue
        
        # Use the first candidate carbon: check if it is attached to a long aliphatic chain.
        # We are only interested if the chain length is at least 12 carbons.
        for c_atom in candidate_c_atoms:
            chain_length = dfs_chain(c_atom, {c_atom.GetIdx()})
            if chain_length >= 12:
                reason = (f"Found sulfonic acid group with S(=O)(=O) and chain length = {chain_length}; "
                          "and polar head (ring) is present")
                return True, reason
            else:
                # Even if a sulfonic acid group is found, if the chain is too short we record that.
                short_reason = f"Found sulfonic acid group, but attached alkyl chain too short (chain length = {chain_length})"
        # If none of the candidate carbons yield a long enough chain for this sulfur,
        # continue searching.
    return False, "No sulfonic acid residue attached to a sufficiently long lipid chain found"