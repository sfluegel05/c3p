"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: Any molecule that contains two amino‐acid residues connected via peptide linkages.
This implementation accepts both linear dipeptides (one backbone peptide bond)
and cyclic dipeptides (diketopiperazines – a 6‐membered ring with 2 peptide bonds).
The strategy is to identify candidate backbone bonds using an amide SMARTS pattern
but then to “break” those bonds and check that exactly two fragments result,
with fragment sizes in a range typical of amino‐acid residues. We also check that 
the overall molecule is small (e.g. heavy atom count below a cutoff). Finally, 
we look for the cyclic dipeptide (diketopiperazine) case as an alternate.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_dipeptide(smiles: str):
    """
    Determines if the given SMILES corresponds to a dipeptide.
    
    For a linear dipeptide, we require that exactly one backbone peptide bond is found.
    This candidate peptide bond is subject to several filters:
      - It must match the basic amide pattern ([C](=O)[N]) and not be inside an aromatic ring.
      - Breaking that bond yields exactly two fragments.
      - Each fragment has a reasonable number of heavy atoms (i.e. between 3 and 40).
    For cyclic dipeptides (diketopiperazines), we look for a 6–membered ring
    that contains exactly two backbone peptide bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a dipeptide, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES and add hydrogens (so that potential alpha-carbons are explicit)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # Check overall size. Typically dipeptides are small molecules.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(heavy_atoms) > 60:
        return False, f"Too many heavy atoms ({len(heavy_atoms)}) for a typical dipeptide."
    
    # Define a basic amide bond pattern without atom mapping (we will extract indices later)
    # We use pattern: carbon with a double-bonded oxygen and a single-bonded nitrogen.
    amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")
    matches = mol.GetSubstructMatches(amide_pattern)
    if not matches:
        return False, "No amide (peptide) bond pattern detected."
    
    valid_linear_bonds = []
    # For each match, extract candidate bond and apply filters
    for match in matches:
        # match is a tuple of atom indices; by our SMARTS the first is carbonyl carbon and last is N.
        c_idx = match[0]
        n_idx = match[-1]
        # Get the bond between these atoms
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond is None:
            continue
        # Skip candidate if either atom is aromatic (often artifacts in non-peptide amides)
        if mol.GetAtomWithIdx(c_idx).GetIsAromatic() or mol.GetAtomWithIdx(n_idx).GetIsAromatic():
            continue
        # For a linear dipeptide we do not want the bond to be in a ring.
        if bond.IsInRing():
            # We'll treat cyclic candidates later.
            continue
        # Now “break” the bond and see if we get exactly two fragments.
        try:
            # FragmentOnBonds will add dummy atoms ([*]) at the break sites.
            frag_mol = rdmolops.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
        except Exception as e:
            continue
        frags = Chem.GetMolFrags(frag_mol, asMols=True)
        if len(frags) != 2:
            continue
        # For each fragment, count heavy atoms (non-hydrogen).
        frag_ok = True
        for frag in frags:
            frag_heavy = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1]
            # Typical amino acid residue (with protecting groups removed) has at least 3 heavy atoms,
            # and in dipeptides residues are usually not huge; set an arbitrary upper bound.
            if not (3 <= len(frag_heavy) <= 40):
                frag_ok = False
                break
        if frag_ok:
            valid_linear_bonds.append((c_idx, n_idx))
    
    # If exactly one valid non‐ring peptide bond candidate is found then we classify as a linear dipeptide.
    if len(valid_linear_bonds) == 1:
        return True, "Dipeptide detected: a single peptide bond connects two amino acid residues."
    elif len(valid_linear_bonds) > 1:
        return False, f"Found {len(valid_linear_bonds)} valid peptide bond candidates; ambiguous for a dipeptide."
    
    # Next, check for cyclic dipeptide (diketopiperazine).
    # In a diketopiperazine, the two peptide bonds appear in a 6–membered ring.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        # Count candidate peptide bonds whose both atoms lie in the ring.
        bond_count_in_ring = 0
        for match in matches:
            c_idx = match[0]
            n_idx = match[-1]
            # Get the bond
            bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
            if bond is None:
                continue
            # For cyclic dipeptide, we expect the two amide bonds to be part of the same 6-membered ring.
            if c_idx in ring and n_idx in ring:
                # Also skip aromatic bonds even in rings.
                if mol.GetAtomWithIdx(c_idx).GetIsAromatic() or mol.GetAtomWithIdx(n_idx).GetIsAromatic():
                    continue
                bond_count_in_ring += 1
        if bond_count_in_ring == 2:
            return True, "Cyclic dipeptide (diketopiperazine) detected: 6-membered ring with 2 peptide bonds."
    
    # If we reach here, no valid linear or cyclic dipeptide candidate is found.
    return False, "No valid peptide bond connecting two appropriate fragments was found."

# Example usage for testing:
if __name__ == "__main__":
    test_smiles_list = [
        # These examples are known dipeptides (linear or cyclic)
        "S(=O)(=N)(CC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)O",  # gamma-glutamylmethionine sulfoximine
        "SCC(NC(=O)C(N)C(O)C)C(O)=O",  # Threoninyl-Cysteine
        "CC(C)[C@H](NC(=O)[C@H](C)N)C(O)=O",  # Ala-Val
        "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)NCC(O)=O",  # Trp-Gly
        "C1C[C@@]2([H])C(N[C@H](C(N2C1)=O)CC(C)C)=O",  # cyclo(L-Leu-L-Pro) cyclic dipeptide example
        # A molecule that is too large or ambiguous:
        "NCC(N1[C@@H](CCC1)C(=O)NC=2C=CC3=C(C2)C=CC=C3)=O"  # glycyl-L-proline 2-naphthylamide (ambiguous)
    ]
    for sm in test_smiles_list:
        res, reason = is_dipeptide(sm)
        print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")