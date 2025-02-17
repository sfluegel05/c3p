"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Defined as: any member of the class of chlorobenzenes carrying three chloro substituents 
at unspecified positions.
The improved classifier now (i) rejects molecules with multiple fragments (i.e. salts),
(ii) rejects molecules with metal atoms, and 
(iii) for each aromatic six‐membered ring it only counts it as a trichlorobenzene
if the ring has exactly three chlorine atoms as substituents and all other substituents
are attached via carbon (atomic number 6) or oxygen (atomic number 8).
"""

from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if the molecule is a trichlorobenzene.
    A trichlorobenzene is defined as a molecule containing a six-membered aromatic ring 
    that has exactly three chlorine atoms attached as substituents and that, aside from
    chlorine, all other substituents attached directly to the ring atoms come from either 
    carbon or oxygen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a trichlorobenzene, False otherwise.
        str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains more than one fragment (e.g. salts or counter ions)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Molecule has multiple fragments (possible salt or counter ion present)"
    # (Alternatively, one could try to work only with the largest fragment.)

    # Reject molecules containing metal atoms.
    for atom in mol.GetAtoms():
        # The IsMetal() function is available on rdkit atom objects.
        if atom.GetIsMetal():
            return False, "Molecule contains metal atoms"
    
    # Get ring information 
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Loop over each ring.
    for ring in atom_rings:
        # We want 6-membered rings.
        if len(ring) != 6:
            continue
        
        # Check if all atoms in the ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        chloro_count = 0
        ring_ok = True  # flag for allowed substituent types on the ring
        
        # For each atom in the ring check its neighbors that are not part of the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider substituents not in the ring.
                if nbr.GetIdx() in ring:
                    continue
                atomic_num = nbr.GetAtomicNum()
                if atomic_num == 17:
                    chloro_count += 1
                else:
                    # If neighbor is not Cl, only allow attachments via carbon or oxygen.
                    if atomic_num not in (6, 8):
                        ring_ok = False
                        # Break out early on this atom since we found a disallowed substituent.
                        break
            if not ring_ok:
                break

        if not ring_ok:
            # This ring has a disallowed substituent type, so do not classify it.
            continue

        # If found exactly three chlorine substituents on this aromatic ring then it qualifies.
        if chloro_count == 3:
            return True, "Found a benzene ring with exactly three chloro substituents (all other substituents attached via C or O)"
    
    return False, "No qualifying aromatic benzene ring with exactly three chloro substituents found"


# Example usage (these prints can be removed or commented out when integrating into a larger system):
if __name__ == "__main__":
    test_smiles = [
        "Clc1cccc(Cl)c1Cl",  # 1,2,3-trichlorobenzene; valid
        "C1(=CC=C(C(=C1Cl)CC([O-])=O)Cl)Cl.[Na+]",  # chlorfenac-sodium; should be rejected
        "CON(C(C)Cc1c(Cl)cc(Cl)cc1Cl)C(=O)c1cn(C)nc1C(F)F"  # valid 3-(difluoromethyl)-...pyrazole carboxamide
    ]
    for smi in test_smiles:
        result, reason = is_trichlorobenzene(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")