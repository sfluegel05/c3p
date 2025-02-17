"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide
Definition: Any cyclic ether in which the oxygen atom forms part of a 3‐membered ring.
That is, a 3‐membered ring (all single bonds) containing exactly one oxygen and two carbons,
where the carbons are sp³ and no exocyclic double bond (i.e. conjugation) is present on the epoxide carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_epoxide(smiles: str):
    """
    Determines whether the molecule (given as a SMILES string) can be classified as an epoxide.
    We define an epoxide as any 3-membered ring containing exactly one oxygen and two carbons,
    where all bonds in the ring are single, the carbons are sp3-hybridized and no exocyclic double bond 
    is attached to the ring carbons.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if at least one valid epoxide ring is found, False otherwise.
        str: Detailed reason for the classification.
    """
    # Try to parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a 3-membered ring where the ring atoms are: carbon - oxygen - carbon.
    # The '-' bonds enforce that the bonds are single.
    # We also require that the atoms belong to a 3-membered ring (r3).
    epoxide_smarts = "[#6;r3]-[#8;r3]-[#6;r3]"
    pattern = Chem.MolFromSmarts(epoxide_smarts)
    if pattern is None:
        return False, "SMARTS pattern error"
    
    # Find all matches (each match is a tuple of three atom indices that match the pattern).
    matches = mol.GetSubstructMatches(pattern, uniquify=True)
    
    valid_count = 0
    reasons = []
    
    # For each match, perform additional sanity checks.
    for match in matches:
        # Get the atoms involved in the potential epoxide ring.
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        # Quick check: if any atom is aromatic, skip this match.
        if any(atom.GetIsAromatic() for atom in atoms):
            continue
        
        # Verify that the atoms indeed form a ring (their three bonds should connect cyclically).
        # (The SMARTS ensures connectivity but we will explicitly check the bonds in our match.)
        # Order the match atoms in a ring order using the pattern bond connectivity.
        # We check bonds between: atom1-atom2, atom2-atom3, and atom3-atom1.
        if not (mol.GetBondBetweenAtoms(match[0], match[1]) and
                mol.GetBondBetweenAtoms(match[1], match[2]) and
                mol.GetBondBetweenAtoms(match[2], match[0])):
            continue
        
        # Check that all bonds in the ring are single (the '-' SMARTS symbol should ensure that,
        # but we revalidate here).
        bond_indices = [(match[0], match[1]), (match[1], match[2]), (match[2], match[0])]
        single_bonds = True
        for i, j in bond_indices:
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                single_bonds = False
                break
        if not single_bonds:
            continue
        
        # Ensure that the ring contains exactly one oxygen and two carbons.
        oxygen_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if oxygen_count != 1 or carbon_count != 2:
            continue
        
        # For each carbon in the ring, check that it is sp3 hybridized.
        # This reduces false positives where the carbon might be vinylic or otherwise conjugated.
        carbons_sp3 = True
        for atom in atoms:
            if atom.GetAtomicNum() == 6:
                if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    carbons_sp3 = False
                    break
        if not carbons_sp3:
            continue
        
        # Check that none of the carbon atoms in the ring is directly attached to an exocyclic double bond.
        # (i.e. any bond from a ring carbon to an atom outside the epoxide ring should not be a double bond.)
        ring_atom_set = set(match)
        exocyclic_double = False
        for atom in atoms:
            if atom.GetAtomicNum() == 6:  # only for carbons in the ring
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetIdx() in ring_atom_set:
                        continue
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        exocyclic_double = True
                        break
                if exocyclic_double:
                    break
        if exocyclic_double:
            continue

        # If all criteria are satisfied then count this epoxide ring.
        valid_count += 1

    if valid_count > 0:
        return True, f"Found epoxide ring(s): {valid_count} occurrence(s) of a 3-membered cyclic ether"
    else:
        return False, "No valid epoxide ring (3-membered cyclic ether with one oxygen and two sp3 carbons and no exocyclic double bonds) found"


# (Optional test code)
if __name__ == "__main__":
    # A few test cases to check the classifier.
    test_smiles = [
        # True positives
        "CCCCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O",     # juvenile hormone I acid
        "Fc1ccc(cc1)[C@]1(Cn2cncn2)O[C@H]1c1ccccc1Cl",          # (2R,3S)-epoxiconazole
        "C1O[C@@H]1c1ccccc1",                                   # (R)-styrene oxide
        # A couple that were problematic
        "O1C(C1C#CC#CC2OC2CO)CO",                               # Repandiol (has 2 epoxide rings)
        "N[C@@H](C[C@@H]1CCC(=O)[C@@H]2O[C@H]12)C(O)=O",        # anticapsin
    ]
    for s in test_smiles:
        result, reason = is_epoxide(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")