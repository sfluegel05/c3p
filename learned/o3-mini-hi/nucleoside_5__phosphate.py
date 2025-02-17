"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate
Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which the
sugar’s C-5 (the exocyclic carbon, typically CH2OH in ribose/deoxyribose) is phosphorylated
(mono-, di-, tri- or tetra-phosphate).

Heuristics in this improved version:
  1. Reject very heavy molecules (MW >800 Da) because nucleosides are of moderate size.
  2. Identify a nucleobase candidate as any fully aromatic ring (all atoms aromatic) having ≥2 N.
  3. Identify a sugar candidate by looking for a 5‐membered ring (furanose) that contains 4 carbons and 1 oxygen.
     • Exactly one atom in the ring must be linked to the nucleobase (the anomeric carbon).
     • One of the other ring carbons (preferably a primary –CH2– group, which has at least 2 hydrogens)
       must have an exocyclic oxygen substituent that in turn is attached to a phosphorus atom.
  4. If no closed sugar candidate is found, a fallback open‐chain pattern [CH2]-[O]-[P] is used,
     provided that the CH2 is attached to the nucleobase.
     
If these conditions are met, the function returns True and a reason; otherwise, it returns False with the reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the nucleoside 5'-phosphate criteria, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the input SMILES and add explicit hydrogens to help count them.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol_with_H = Chem.AddHs(mol)
    
    # Reject overly heavy molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} is too high for a typical nucleoside phosphate (<800 Da)."
    
    # Identify nucleobase candidate atoms.
    # A nucleobase candidate is defined as any fully aromatic ring with at least 2 nitrogen atoms.
    nucleobase_atoms = set()
    ring_info = mol_with_H.GetRingInfo()
    for ring in ring_info.AtomRings():
        if all(mol_with_H.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            n_nitrogens = sum(1 for i in ring if mol_with_H.GetAtomWithIdx(i).GetAtomicNum() == 7)
            if n_nitrogens >= 2:
                nucleobase_atoms.update(ring)
    if not nucleobase_atoms:
        return False, "No nucleobase candidate found (aromatic ring with ≥2 nitrogen atoms)."
    
    # Search for a closed sugar candidate:
    # Look for a 5-membered ring with 4 carbons and 1 oxygen.
    # It must connect to the nucleobase in exactly one position (the anomeric carbon).
    sugar_ring_found = None
    anomeric_idx = None
    phosphate_found = False
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        oxygens = [i for i in ring if mol_with_H.GetAtomWithIdx(i).GetAtomicNum() == 8]
        carbons = [i for i in ring if mol_with_H.GetAtomWithIdx(i).GetAtomicNum() == 6]
        if len(oxygens) != 1 or len(carbons) != 4:
            continue
        
        # Determine which atom(s) in this ring attach to a nucleobase candidate.
        base_connections = []
        for idx in ring:
            atom = mol_with_H.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # use idx from mol_with_H (neighbors are from the same molecule)
                if nbr.GetIdx() in nucleobase_atoms:
                    base_connections.append(idx)
                    break
        if len(base_connections) != 1:
            continue  # avoid rings attached more than one time (possible dinucleotides)
        # The single connection is the anomeric carbon.
        anomeric_idx = base_connections[0]
        
        # Now check for an exocyclic phosphoryl group on the sugar.
        # Look among the non-anomeric ring atoms for one that resembles a CH2 (at least 2 H)
        # and which has a neighbor oxygen that is bound to a phosphorus.
        for idx in ring:
            if idx == anomeric_idx:
                continue
            atom = mol_with_H.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() >= 2:  # candidate for primary CH2
                for nbr in atom.GetNeighbors():
                    # Ensure this neighbor is exocyclic (i.e. not in the sugar ring)
                    if nbr.GetIdx() in ring:
                        continue
                    if nbr.GetAtomicNum() == 8:  # potential oxygen substituent
                        # Check if this oxygen is bound to a phosphorus
                        for o_nbr in nbr.GetNeighbors():
                            if o_nbr.GetAtomicNum() == 15:
                                phosphate_found = True
                                sugar_ring_found = set(ring)
                                break
                    if phosphate_found:
                        break
            if phosphate_found:
                break
        if sugar_ring_found and phosphate_found:
            break
    
    # Fallback: if no closed sugar ring was found, search for an open‐chain fragment.
    if not sugar_ring_found or not phosphate_found:
        # Look for a fragment matching [CH2]-[O]-[P]
        patt = Chem.MolFromSmarts("[CH2]-[O]-[P]")
        matches = mol_with_H.GetSubstructMatches(patt)
        for match in matches:
            c_idx, o_idx, p_idx = match
            # Check that the CH2 atom is attached to a nucleobase.
            for nbr in mol_with_H.GetAtomWithIdx(c_idx).GetNeighbors():
                if nbr.GetIdx() in nucleobase_atoms:
                    phosphate_found = True
                    sugar_ring_found = {c_idx}  # mark this as our sugar candidate fragment
                    break
            if phosphate_found:
                break

    if not sugar_ring_found or not phosphate_found:
        return False, ("No sugar candidate found that is attached to a nucleobase and bears an exocyclic "
                       "oxygen bound to phosphorus at the C-5 position.")
    
    # Finally, check that the sugar candidate does not have extra connections to the nucleobase.
    # For a bona fide nucleoside, the sugar should attach to the nucleobase in exactly one place.
    connection_count = 0
    for idx in sugar_ring_found:
        for nbr in mol_with_H.GetAtomWithIdx(idx).GetNeighbors():
            if nbr.GetIdx() in nucleobase_atoms:
                connection_count += 1
    if connection_count != 1:
        return False, ("Sugar candidate has an incorrect number of nucleobase connections ("
                       f"{connection_count} found; expected exactly one), possibly indicating a dinucleotide or over‐phosphorylated species.")
    
    return True, ("Molecule contains a nucleobase (aromatic heterocycle with ≥2 nitrogen atoms) attached to "
                  "a sugar candidate (a 5-membered furanose or open-chain CH2 fragment) that bears a phosphate group "
                  "at the likely 5'-position.")

# (Optional) Testing examples:
# test_smiles = [
#     "CC(=O)Nc1ccn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)n1",  # N(4)-acetylcytidine 5'-monophosphate (expected True)
#     "C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O)",  # FMN-L-threonine (True)
#     "NC1=NC(=O)N(C=C1)[C@H]1C[C@H](O)[C@@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)O1",  # dCTP (should be borderline True/False)
# ]
# for smi in test_smiles:
#     res, msg = is_nucleoside_5__phosphate(smi)
#     print(smi, "==>", res, msg)