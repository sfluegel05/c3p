"""
Classifies: CHEBI:33567 catecholamine
"""
#!/usr/bin/env python
"""
Classifies: catecholamine (4-(2-Aminoethyl)pyrocatechol and its derivatives)
A catecholamine is defined by the presence of a catechol moiety (benzene ring with two adjacent hydroxyl groups)
and an aminoethyl side chain (i.e. –CH2–CH2–NH2 or a substituted variant attached to the catechol ring).
This program uses heuristic substructure detection to classify a SMILES string.
"""

from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    The heuristic rules are:
      1. The molecule must contain at least one benzene ring (6-membered aromatic ring)
         that has two adjacent ring atoms substituted with hydroxyl (-OH) groups.
      2. The molecule must contain a two-carbon (aminoethyl) side chain ending with a nitrogen 
         (allowing for common substitutions) attached to the aromatic system.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a catecholamine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flag to mark if a catechol (benzene-1,2-diol) moiety has been found.
    catechol_found = False
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    
    # Iterate over all rings in the molecule.
    for ring in ring_info.AtomRings():
        # Look only at rings that are 6-membered.
        if len(ring) == 6:
            # Check that every atom in the ring is aromatic.
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            
            # For atoms in the ring, collect indices that have a hydroxyl (-OH) group.
            oh_atoms = []
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Check neighbors: if any neighbor is oxygen, not in the same ring,
                # and that oxygen is attached only to this atom then it is likely an OH.
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in ring):
                        # Heuristic: an oxygen with a single heavy-atom neighbor is likely -OH.
                        if nbr.GetDegree() == 1:
                            oh_atoms.append(idx)
                            break
            # Now check if there are at least two adjacent ring atoms (cyclically adjacent)
            # that carry an OH group.
            if len(oh_atoms) >= 2:
                ring_length = len(ring)
                for i in range(ring_length):
                    # Consider adjacent atoms in the ring (cyclically adjacent)
                    current_atom = ring[i]
                    next_atom = ring[(i+1) % ring_length]
                    if current_atom in oh_atoms and next_atom in oh_atoms:
                        catechol_found = True
                        break
            if catechol_found:
                break  # Stop if we have found a valid catechol moiety.
    
    if not catechol_found:
        return False, "No catechol moiety (benzene ring with adjacent hydroxyl groups) found"
    
    # Check for an aminoethyl side chain.
    # SMARTS pattern: "cCCN" represents an aromatic carbon connected to two aliphatic carbons ending in a nitrogen.
    aminoethyl_pattern = Chem.MolFromSmarts("cCCN")
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl side chain (pattern 'cCCN') found"
    
    return True, "Contains a catechol moiety with an aminoethyl side chain"
    
# Example test cases (uncomment to run tests):
# test_smiles = [
#     "NCCc1ccc(O)c(O)c1",  # dopamine
#     "C[C@H](N)[C@H](O)c1ccc(O)c(O)c1",  # (-)-alpha-Methylnoradrenaline
#     "CC1=C(C(=CC(=C1O)O)[N+](=O)[O-])CCN"  # 4-(2-aminoethyl)-5-nitrobenzene-1,2-diol
# ]
# for smi in test_smiles:
#     result, reason = is_catecholamine(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")