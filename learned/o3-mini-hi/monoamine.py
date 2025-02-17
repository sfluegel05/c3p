"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: Monoamine, defined as an aralkylamino compound which contains one amino group 
connected to an aromatic ring by a two‐carbon chain.
Examples include tyramine, (R)-noradrenaline, dopamine, etc.
"""

from rdkit import Chem
from rdkit.Chem.rdchem import HybridizationType

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an aralkylamino compound having one amino group connected
    to an aromatic ring by a two‐carbon chain. In this implementation we search for a nitrogen (not aromatic)
    that is not involved in an amide-like bond and which has a shortest path of exactly 3 bonds to an
    aromatic carbon. The two intervening atoms must be non‐aromatic carbons (and ideally sp3) forming the
    two-carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a monoamine, False otherwise.
        str: A message detailing the reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    valid_monoamine_count = 0

    # Iterate over every atom. We look for non-aromatic nitrogen atoms that are not part of an amide.
    for atom in mol.GetAtoms():
        # Look at nitrogen atoms (atomic number 7), which are not aromatic.
        if atom.GetAtomicNum() != 7 or atom.GetIsAromatic():
            continue

        # Exclude nitrogen that might be in an amide (or similar) environment:
        # If a neighbor (typically a carbon) is double bonded to an oxygen, then skip.
        is_amide = False
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 6:  # carbon neighbor
                for bond in nb.GetBonds():
                    # Check for a double bond to an oxygen.
                    if bond.GetBondTypeAsDouble() == 2:
                        other = bond.GetOtherAtom(nb)
                        if other.GetAtomicNum() == 8:
                            is_amide = True
                            break
                if is_amide:
                    break
        if is_amide:
            continue

        # For this nitrogen, search for an aromatic carbon reachable by a path of exactly 3 bonds.
        valid_substruct = False
        for target in mol.GetAtoms():
            # Consider only aromatic carbons.
            if target.GetAtomicNum() != 6 or not target.GetIsAromatic():
                continue
            # Avoid trivial self-match.
            if target.GetIdx() == atom.GetIdx():
                continue

            # Compute the shortest bond path between the nitrogen and the target aromatic carbon.
            path = Chem.rdmolops.GetShortestPath(mol, atom.GetIdx(), target.GetIdx())
            # We require exactly 4 atoms in the path (i.e. 3 bonds).
            if len(path) != 4:
                continue

            # Check that the two atoms in between are carbons, are non‐aromatic, and are SP3.
            intermediate1 = mol.GetAtomWithIdx(path[1])
            intermediate2 = mol.GetAtomWithIdx(path[2])
            if (intermediate1.GetAtomicNum() == 6 and intermediate2.GetAtomicNum() == 6 and
                not intermediate1.GetIsAromatic() and not intermediate2.GetIsAromatic() and
                intermediate1.GetHybridization() == HybridizationType.SP3 and
                intermediate2.GetHybridization() == HybridizationType.SP3):
                valid_substruct = True
                break

        # If we found the expected connectivity on this nitrogen, count it.
        if valid_substruct:
            valid_monoamine_count += 1

    # The molecule should contain exactly one valid monoamine group.
    if valid_monoamine_count == 0:
        return False, "No aralkylamino moiety (amine connected by a two‐carbon chain to an aromatic ring) found"
    elif valid_monoamine_count > 1:
        return False, f"Found {valid_monoamine_count} aralkylamino moieties; expected exactly one monoamine group"
    else:
        return True, "Contains exactly one amino group connected via a two‐carbon chain to an aromatic ring"

# Example usage (for manual testing)
# test_smiles = [
#   "CNC[C@@H](O)c1ccc(O)c(O)c1",  # (S)-adrenaline
#   "NCCc1ccc(O)cc1",              # tyramine
#   "[NH3+]CCc1ccc(O)cc1",          # tyraminium (should be accepted)
#   "OC(=O)[C@H](Cc1ccc(O)c(O)c1)\\N=C/C=C1C[C@H](NC(=C\\1)C(O)=O)C(O)=O"  # Dopaxanthin (likely to be rejected)
# ]
# for s in test_smiles:
#   result, reason = is_monoamine(s)
#   print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")