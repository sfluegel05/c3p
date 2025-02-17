"""
Classifies: CHEBI:48039 dihydroflavonols
"""
#!/usr/bin/env python3
"""
Classifies: dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.
This implementation looks for the key flavanone (2-phenylchroman-4-one) scaffold by identifying:
  - a carbonyl group (C=O) directly bonded to a saturated carbon that bears an –OH (expected to be at C3),
  - verifying that these atoms belong to a six-membered ring that also contains a heterocyclic oxygen (expected at C1),
  - and ensuring that a different ring carbon (expected C2) is bonded to an external aromatic (B‐ring) substituent.
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a given SMILES string corresponds to a dihydroflavonol.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule conforms to the dihydroflavonol (2-phenylchroman-4-one with C3-OH) scaffold.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, search for the key fragment: Carbonyl (C=O) attached to a saturated carbon with an OH.
    # We relax the pattern by using "[CH]" instead of "[C]" so that chiral annotations (like [C@@H]) are tolerated.
    query = Chem.MolFromSmarts("C(=O)[CH]([OH])")
    if query is None:
        return None, None  # Should not occur.
    
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Molecule does not contain the key fragment C(=O)[CH]([OH]) expected in dihydroflavonols (indicative of C4 and C3)"
    
    # Get all rings (as sets of atom indices) found by RDKit.
    rings = [set(r) for r in Chem.GetSymmSSSR(mol)]
    
    # Loop over each fragment match.
    for match in matches:
        # In our query:
        #   match[0] is the carbonyl carbon (C4)
        #   match[1] is the saturated carbon bearing -OH (C3)
        carbonyl_idx = match[0]
        c3_idx = match[1]
        
        # Find six-membered rings that contain both the carbonyl and C3 atoms.
        candidate_rings = [ring for ring in rings if (carbonyl_idx in ring and c3_idx in ring and len(ring) == 6)]
        if not candidate_rings:
            continue  # Try the next fragment match if this one is not in any six-membered ring.
        
        # Evaluate each candidate ring.
        for ring in candidate_rings:
            # (i) Verify that the ring contains a heterocyclic oxygen, expected to be part of the chroman ring.
            # We only consider ring oxygens that are likely not -OH (i.e. have more than one connection).
            ring_oxygens = []
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8 and atom.GetSymbol() == "O" and atom.GetDegree() > 1:
                    ring_oxygens.append(idx)
            if not ring_oxygens:
                continue  # No suitable ring oxygen found in this candidate ring.
            
            # (ii) Look for a ring carbon (excluding the carbonyl and C3) that is bonded externally to an aromatic system.
            aromatic_attachment_found = False
            for idx in ring:
                if idx in (carbonyl_idx, c3_idx):
                    continue  # Skip atoms already used in the key fragment.
                atom = mol.GetAtomWithIdx(idx)
                # Consider only carbon atoms.
                if atom.GetAtomicNum() != 6:
                    continue
                # Optionally, one might check if this atom is adjacent to the ring oxygen.
                # Here we require that at least one neighbor (outside the ring) is aromatic.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ring and nbr.GetIsAromatic():
                        aromatic_attachment_found = True
                        break
                if aromatic_attachment_found:
                    break  # Found an atom meeting the condition.
            
            if aromatic_attachment_found:
                return True, ("Matches dihydroflavonol core: contains a 2-phenylchroman-4-one scaffold with an -OH at the C3 position "
                              "and a ring oxygen in a six-membered heterocycle.")
    
    return False, "Molecule does not match the dihydroflavonol core criteria."

# Example usage (for testing purposes)
if __name__ == "__main__":
    # List of some example SMILES strings for dihydroflavonols.
    examples = [
        "O1C(C(O)C(=O)C=2C1=CC=3OCOC3C2O",  # 3,5-Dihydroxy-6,7-methylenedioxyflavanone
        "O[C@@H]1[C@H](Oc2cc(O)ccc2C1=O)c1ccc(O)cc1",    # garbanzol
        "CC(C)=CCc1c(O)cc(O)c2C(=O)[C@H](O)[C@H](Oc12)c1ccccc1",  # glepidotin B
        "O1[C@H]([C@H](O)C(=O)C=2C1=CC(O)=CC2O)C3=CC=C(O)C=C3",  # (2S,3S)-3,5,7-Trihydroxy-2-(4-hydroxyphenyl)-2,3-dihydro-4H-chromen-4-one
        "COC1=C(O)C=CC(=C1)[C@H]1OC2=CC(O)=CC(O)=C2C(=O)[C@@H]1O", # dihydroisorhamnetin
        "OC1C(Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",           # (+)-taxifolin
        "O[C@H]1[C@@H](Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1",  # (-)-dihydromyricetin
    ]
    for s in examples:
        result, reason = is_dihydroflavonols(s)
        print(f"SMILES: {s}\nClassified as dihydroflavonol? {result}\nReason: {reason}\n")