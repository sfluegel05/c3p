"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl derivatives
Definition: A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives.
A true para-terphenyl should contain a central benzene ring (six-membered, aromatic, only carbons)
with two benzene rings attached via single (non-fused) bonds at the para positions (i.e. atoms separated by 3 positions).
"""

from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl derivative based on its SMILES string.

    The algorithm works as follows:
      1. Parse the SMILES with RDKit and sanitize the molecule.
      2. Compute the SSSR (smallest set of smallest rings) and select six-membered rings.
      3. Identify candidate central rings: those that are aromatic, six-membered and composed entirely of carbons.
      4. For each candidate, use the ring ordering provided by SSSR and check each atom for external substituents.
         A valid substituent is a benzene ring (six-membered aromatic, all-carbon) that is attached via a single bond,
         meaning the external ring shares exactly one atom with the central ring.
      5. If at least two substituents are found on the central ring, check if any two are para relative to each other
         (i.e. their positions in the ring differ by 3 modulo 6).
      6. If such a pair is found, the molecule is classified as a para-terphenyl derivative.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule displays a para-terphenyl motif, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error sanitizing molecule: {e}"

    # Use RDKit's SSSR to get ring information (each ring is in cyclic order)
    sssr = Chem.GetSymmSSSR(mol)
    # Filter to six-membered rings; each ring is represented as a tuple of atom indices.
    six_membered_rings = [list(r) for r in sssr if len(r) == 6]
    if not six_membered_rings:
        return False, "No six-membered rings found in the molecule"

    # Identify candidate central rings:
    # Must be aromatic and contain only carbon atoms.
    candidate_central = []
    for ring in six_membered_rings:
        is_candidate = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() and atom.GetSymbol() == "C"):
                is_candidate = False
                break
        if is_candidate:
            candidate_central.append(ring)
    if not candidate_central:
        return False, "No pure carbon aromatic six-membered rings found as candidates for central rings"

    # Now check each candidate central ring for two external benzene substituents in para positions.
    # Among the six-membered rings, external benzene rings are similarly defined.
    valid_external_rings = []
    for ring in six_membered_rings:
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() and atom.GetSymbol() == "C"):
                is_benzene = False
                break
        if is_benzene:
            valid_external_rings.append(set(ring))  # store as a set for overlap tests

    # Loop over candidate central rings.
    for central_ring in candidate_central:
        # The order of atoms in central_ring is given by the SSSR; we assume this ordering is cyclic.
        substituent_positions = []
        n = 6  # six atoms in a benzene
        # For each atom in the central ring, look for an external attachment that qualifies.
        for pos, atom_idx in enumerate(central_ring):
            central_atom = mol.GetAtomWithIdx(atom_idx)
            # Check each neighbor of central_atom
            for nbr in central_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in central_ring:
                    continue  # skip atoms belonging to the central ring
                # Now check: does this neighbor belong to any valid benzene ring that is attached in a non-fused manner?
                for ext_ring in valid_external_rings:
                    # The external ring must share exactly one atom with the central ring (the connection point).
                    if nbr_idx in ext_ring and len(set(central_ring).intersection(ext_ring)) == 1:
                        substituent_positions.append(pos)
                        # Once one substituent is found at this position, no need to search further.
                        break
                else:
                    continue
                break
                    
        # Need at least two substituents for a para-terphenyl skeleton.
        if len(substituent_positions) < 2:
            continue
        # Check if any substituent pair are para on the six-membered ring (positions differing by 3 modulo 6).
        for i in range(len(substituent_positions)):
            for j in range(i+1, len(substituent_positions)):
                diff = abs(substituent_positions[i] - substituent_positions[j])
                circular_diff = min(diff, n - diff)
                if circular_diff == 3:
                    return True, ("Found para-terphenyl skeleton: central benzene ring with two substituent benzene rings "
                                  "attached at para positions (separated by 3 atoms in the ring) via single (non-fused) bonds.")
    return False, "Could not find a central benzene ring with two para-attached benzene substituents (non-fused)."

# Example usage:
# Uncomment the following lines to test with some of the provided SMILES examples:
# test_smiles = [
#     "O=C(OC1=C(OC(=O)C)C(=C(O)C(=C1C2=CC(O)=C(O)C=C2)O)C3=CC=C(O)C=C3)C",  # 2',3'-diacetoxy-3,4,5',6',4''-pentahydroxy-p-terphenyl
#     "COc1cc(-c2ccccc2)c(OC)c(O)c1-c1ccc(O)c(O)c1",  # A simplified p-terphenyl derivative
# ]
# for s in test_smiles:
#     result, reason = is_para_terphenyl(s)
#     print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")