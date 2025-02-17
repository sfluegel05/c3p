"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: Secondary amine 
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
Only exocyclic secondary amine groups (i.e. the candidate nitrogen must not be part of a ring)
are counted. In addition, if one of the “hydrogen positions” is replaced by a nitroso group (-N=O),
that replacement is allowed.

This implementation no longer uses fixed SMARTS but instead iterates over each nitrogen atom 
in the molecule and checks:
  1. The nitrogen is not in a ring.
  2. It is not directly bound to a carbonyl group (i.e. no neighboring carbon has a double-bonded oxygen).
  3. It has two hydrocarbyl substituents. In our case we require that a substituent is a carbon atom.
  4. The remaining hydrogen count is exactly 1 (normal secondary amine) or 0 when one substituent 
     is a nitroso group. A nitroso group is detected when one neighboring atom is itself a nitrogen 
     that is double-bonded to an oxygen.
     
If a candidate nitrogen meets either criteria, the function returns True plus an explanation.
Otherwise it returns False and an explanation.
"""

from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if the given molecule (SMILES string) contains an exocyclic secondary amine group.
    A secondary amine herein is defined as a nitrogen derived from ammonia by replacing exactly two hydrogens 
    with hydrocarbyl groups. In one variant, one of these hydrogens may be replaced by a nitroso (-N=O) group.
    
    The candidate nitrogen must not be part of any ring and must not be directly bonded to a carbonyl carbon.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if a qualifying secondary amine is found; False otherwise.
        str: Explanation message.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms; we are interested in nitrogen (atomic number 7)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        
        # Exclude nitrogen atoms that are part of any ring
        if atom.IsInRing():
            continue
        
        # Also, for our purposes the nitrogen should not be aromatic.
        # (Most exocyclic amines are drawn non-aromatically.)
        if atom.GetIsAromatic():
            continue

        # Get total number of attached hydrogens (implicit+explicit)
        n_H = atom.GetTotalNumHs()
        
        # Initialize counters for substituents:
        # carbon_count counts substituents that are hydrocarbyl (here, direct neighbor must be carbon)
        # nitroso_count counts substituents that qualify as a nitroso group.
        carbon_count = 0
        nitroso_count = 0
        
        # Iterate over the neighbors of the candidate nitrogen.
        for nbr in atom.GetNeighbors():
            # If the neighbor is a carbon, count it as a hydrocarbyl substituent.
            if nbr.GetAtomicNum() == 6:
                carbon_count += 1
            # Check if the neighbor could be a nitroso substituent.
            # A nitroso group will be a nitrogen that is double-bonded to an oxygen.
            elif nbr.GetAtomicNum() == 7:
                found_double_O = False
                for bond in nbr.GetBonds():
                    other = bond.GetOtherAtom(nbr)
                    # Check if the neighbor nitrogen is double-bonded to oxygen (atomic num 8).
                    if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() > 1.9:
                        found_double_O = True
                        break
                if found_double_O:
                    nitroso_count += 1
        
        # Exclude candidate nitrogen if it is directly bonded to a carbonyl group.
        # That means: Any carbon neighbor that is double-bonded to an oxygen (C=O) will disqualify.
        carbonyl_attached = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                for bond in nbr.GetBonds():
                    # Check bond type and that the other partner is oxygen.
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() > 1.9:
                        carbonyl_attached = True
                        break
            if carbonyl_attached:
                break
        if carbonyl_attached:
            continue
        
        # Now apply our criteria:
        # (A) A "normal" secondary amine: one hydrogen attached and two carbon substituents.
        if n_H == 1 and carbon_count == 2 and nitroso_count == 0:
            return True, ("Found an exocyclic secondary amine substructure (aliphatic variant): "
                          "a nitrogen with one hydrogen and two hydrocarbyl (carbon) substituents that is not directly "
                          "bonded to a carbonyl group and not part of a ring.")
        # (B) A nitrosated secondary amine: no hydrogen (0) but two carbon substituents and one nitroso substituent.
        if n_H == 0 and carbon_count == 2 and nitroso_count == 1:
            return True, ("Found an exocyclic secondary amine substructure (nitrosated variant): "
                          "a nitrogen with two hydrocarbyl substituents and a nitroso group (replacing a hydrogen) that is not "
                          "directly bonded to a carbonyl group and not part of a ring.")
    
    return False, ("No exocyclic secondary amine detected meeting the criteria: no nitrogen atom outside a ring was found "
                   "with exactly two hydrocarbyl (carbon) substituents and one effective hydrogen (or a nitroso group in place of one hydrogen).")

# For debugging and testing the function with a number of examples:
if __name__ == "__main__":
    test_molecules = [
        ("(R)-dobutamine", "C=1(C(=CC=C(C1)CCN[C@@H](CCC=2C=CC(=CC2)O)C)O)O"),
        ("N-methylcyclohexylamine", "CNC1CCCCC1"),
        ("N(1)-isopropyl-2-methylpropan-1,2-diamine", "CC(C)NCC(C)(C)N"),
        ("(S)-dobutamine", "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1"),
        ("N-(3-aminopropyl)-4-aminobutanal", "[H]C(=O)CCCNCCCN"),
        ("bisoprolol", "CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1"),
        ("N-methylaniline", "CNc1ccccc1"),
        ("nortriptyline", "CNCCC=C1c2ccccc2CCc2ccccc12"),
        ("Bagrelactone A", "O=C1OC2=CC=C([C@H](O)CNC=3C=C1C=CC3O)C=C2"),
        ("methamphetamine", "CN[C@@H](C)Cc1ccccc1"),
        ("dimethylamine", "[H]N(C)C"),
        ("1H-pyrrole", "c1cc[nH]c1")
    ]
    
    for name, smi in test_molecules:
        result, explanation = is_secondary_amine(smi)
        print(f"SMILES: {smi}\n  NAME: {name}\n  -> {result}, Reason: {explanation}\n")