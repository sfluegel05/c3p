"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: Sulfolipid
A sulfolipid is defined as a compound that contains a sulfonic acid residue 
(defined by either S(=O)(=O)[O-] or S(=O)(=O)O) that is attached (directly or via a single oxygen)
to a long, essentially linear, aliphatic chain. Additional criteria:
  - Molecular weight must be above 300 Da.
  - One of the substituents (directly attached or via one bridging oxygen) to the sulfonate 
    must be a sp3 carbon that supports a continuous, mostly unbranched chain of at least 12 sp3 carbons.
If these conditions are met, the molecule is classified as a sulfolipid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    
    A sulfolipid must have a molecular weight above 300 Da, contain a sulfonic acid (or sulfonate) 
    group, and have that sulfonate attached (either directly via a carbon or via one bridging oxygen) 
    to an unbranched aliphatic chain of at least 12 sp3 carbons (counted in a linear walk).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a sulfolipid, False otherwise.
        str: Reason for the classification result.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight (lipid-like molecules are typically large)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for a sulfolipid (mol wt = {mol_wt:.1f})"

    # Define sulfonate SMARTS patterns
    sulfonate_smarts_1 = Chem.MolFromSmarts("S(=O)(=O)[O-]")  # charged sulfonate
    sulfonate_smarts_2 = Chem.MolFromSmarts("S(=O)(=O)O")       # neutral sulfonic acid

    sulfonate_atoms = set()
    for smarts in [sulfonate_smarts_1, sulfonate_smarts_2]:
        matches = mol.GetSubstructMatches(smarts, useChirality=True)
        # In our SMARTS the sulfur is the first atom
        for match in matches:
            sulfonate_atoms.add(match[0])

    if not sulfonate_atoms:
        return False, "No sulfonic acid group detected"

    # Helper: starting from a candidate carbon, walk the chain in one linear direction.
    # We are not following branches – just a simple walk (choosing the longest available path).
    def get_linear_chain_length(atom, came_from=None):
        # Only count if the atom is carbon, sp3, non-aromatic and not in a ring.
        if atom.GetAtomicNum() != 6:
            return 0
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or atom.GetIsAromatic() or atom.IsInRing():
            return 0
        max_length = 1  # count current atom
        # Look for neighbors that are carbons and eligible
        best_extension = 0
        for nbr in atom.GetNeighbors():
            # Do not go back to the atom we came from
            if came_from is not None and nbr.GetIdx() == came_from.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                # Only continue if the neighbor also is sp3, non-aromatic, and acyclic.
                if nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and not nbr.GetIsAromatic() and not nbr.IsInRing():
                    extension = get_linear_chain_length(nbr, atom)
                    if extension > best_extension:
                        best_extension = extension
        return max_length + best_extension

    # For each sulfonate group, check all candidate carbons attached either directly or via a single oxygen.
    for s_idx in sulfonate_atoms:
        sulfur_atom = mol.GetAtomWithIdx(s_idx)
        candidate_c_atoms = []
        # First, check atoms directly bonded to sulfur.
        for nbr in sulfur_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                candidate_c_atoms.append((nbr, sulfur_atom))  # store candidate along with "connector" (here, S)
            # Then, for a bridging oxygen (bond order 1)
            elif nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(sulfur_atom.GetIdx(), nbr.GetIdx())
                if bond is not None and abs(bond.GetBondTypeAsDouble() - 1.0) < 0.01:
                    # Look for a carbon attached to this oxygen (other than the sulfur)
                    for nn in nbr.GetNeighbors():
                        if nn.GetIdx() == sulfur_atom.GetIdx():
                            continue
                        if nn.GetAtomicNum() == 6:
                            candidate_c_atoms.append((nn, nbr))  # candidate reached via oxygen

        # Evaluate each candidate: the candidate carbon should be sp3 and not in a ring (our linear chain walker requires this)
        for cand, connector in candidate_c_atoms:
            if cand.IsInRing() or cand.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # To ensure we count only the aliphatic (lipid‐like) chain, we attempt a linear walk away from cand.
            # If the candidate was reached via a bridging atom, treat that bond as the connection; so we do not walk back.
            chain_length = 0
            # For each eligible neighbor of the candidate (excluding the connector), compute the linear chain length.
            for nbr in cand.GetNeighbors():
                if nbr.GetIdx() == connector.GetIdx():
                    continue
                # Walk the chain starting from nbr (the candidate itself is counted in the chain).
                extension = get_linear_chain_length(nbr, came_from=cand)
                if extension > chain_length:
                    chain_length = extension
            # Add one for the candidate carbon itself.
            total_chain = 1 + chain_length
            if total_chain >= 12:
                reason = (f"Found sulfonic acid group attached to an aliphatic chain (linear chain length = {total_chain}) "
                          f"via atom idx {cand.GetIdx()}")
                return True, reason

    # If none of the candidate carbons yields a linear chain of sufficient length.
    return False, "No sulfonic acid residue found attached to a sufficiently long aliphatic chain"

# Example usage if running this script standalone (for testing):
if __name__ == '__main__':
    test_smiles = [
        # psychosine sulfate (should be sulfolipid)
        "CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO[C@@H]1O[C@H](COS(O)(=O)=O)[C@H](O)[C@H](O)[C@H]1O",
        # A generic aromatic sulfonate (not a sulfolipid)
        "CC1=CC=C(C=C1)S(=O)(=O)O",
        # Example with bridging oxygen (one of the examples provided)
        "S(OC1[C@@H](O)[C@H](O[C@@H](OC[C@H](NC(=O)CCCCCCCCCCCCCCCCCCC)[C@H](O)/C=C/CCCCCCCCCCCCC)C1O)CO)(O)(=O)=O"
    ]
    for smi in test_smiles:
        is_sulfo, reason = is_sulfolipid(smi)
        print(f"SMILES: {smi}\nClassified as sulfolipid? {is_sulfo}\nReason: {reason}\n")