"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary amine
A primary amine is a compound formally derived from ammonia (NH3) by replacing exactly one hydrogen with a hydrocarbyl group.
Thus the defining functional group is –NH2 (with exactly one non‐hydrogen substituent) and not connected to a carbonyl (i.e. not an amide).
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule qualifies as a primary amine based on its SMILES string.
    
    The algorithm works as follows:
      1. Parse the SMILES string and add explicit hydrogens.
      2. Count amide bonds in the molecule (substructure: C(=O)N). If the molecule is both large and polyamide,
         we assume it is a complex multifunctional compound (e.g. a peptide) rather than a simple primary amine.
      3. Iterate over each nitrogen atom (atomic number 7) and, using the explicit hydrogen count,
         check if it is a potential primary amine: exactly two hydrogens and exactly one heavy-atom (non-H) neighbor.
      4. For the candidate nitrogen, check that its single heavy neighbor is not part of a carbonyl (i.e. a C double-bonded to O)
         because that would suggest an amide.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is considered to have a primary amine group, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (ensuring hydrogens are not just implicit).
    mol = Chem.AddHs(mol)
    
    # Count amide bonds using SMARTS. (Pattern: carbonyl directly bonded to a nitrogen.)
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    num_amide = len(amide_matches)
    
    primary_count = 0  # Count candidate primary amine groups
    
    # Iterate over all nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        
        # Get the number of explicit hydrogens on this nitrogen.
        n_explicit_H = atom.GetNumExplicitHs()
        # For a primary amine (R–NH2), we expect exactly two attached hydrogens.
        if n_explicit_H != 2:
            continue
        
        # Get heavy atom (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        
        # The nitrogen in a primary amine should have exactly one non-hydrogen substituent.
        if len(heavy_neighbors) != 1:
            continue
        
        neighbor = heavy_neighbors[0]
        # Exclude if the heavy neighbor is part of a carbonyl:
        # Check if neighbor is double-bonded to an oxygen atom.
        is_adjacent_carbonyl = False
        for bond in neighbor.GetBonds():
            # Check if bond is a double bond.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(neighbor)
                if other.GetAtomicNum() == 8:  # oxygen
                    is_adjacent_carbonyl = True
                    break
        if is_adjacent_carbonyl:
            continue
        
        # If passes all checks, we have found a candidate primary amine group.
        primary_count += 1

    # If no candidate primary amine group was found, report failure.
    if primary_count == 0:
        return False, "No primary amine group (R–NH2) found"
    
    # Apply a heuristic: if the molecule is large and contains multiple amide bonds,
    # its free primary amine may be incidental (e.g., in a peptide).
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_amide > 1 and num_heavy_atoms > 50:
        return False, ("Molecule appears large and contains multiple amide bonds, "
                       "suggesting it is a multifunctional compound rather than a primary amine")
    
    return True, "Contains a primary amine group (R–NH2) defined by exactly one non‐hydrogen substituent and two hydrogens"


# Example usage (for testing a few cases):
if __name__ == "__main__":
    test_smiles = [
        "Nc1ccc(cc1)\\N=N\\c1ccccc1",              # 4-(phenylazo)aniline
        "C[C@@H](N)Cc1ccccc1",                      # (R)-amphetamine
        "C1=CC=CC2=C1C(C3=C(C2=O)C(=C(C(=C3N)O)S(O)(=O)=O)O)=O",  # nuclear fast red free acid
        "CC(C)NCC(C)(C)N",                         # N(1)-isopropyl-2-methylpropan-1,2-diamine
        "CCCCCCCCCCCCCCCCCCN",                      # octadecan-1-amine
        "[H]C(=O)CCCNCCCN",                         # N-(3-aminopropyl)-4-aminobutanal
        "CC(C)C(N)C(C)C",                          # 2,4-dimethylpentan-3-amine
        "CCCCN",                                   # butan-1-amine
        "C1=CC(=C(C=C1O)N)CCN",                     # 3-amino-4-(2-aminoethyl)phenol
        "Nc1ccc2ccccc2c1",                         # 2-naphthylamine
        "CC(C)(N)Cc1ccccc1",                        # phentermine
        "CN",                                      # methylamine
        "COC1=C(O)C=C(CCN)C=C1",                    # 4-methoxytyramine
        "CC(C)(C)NCC(O)c1cc(Cl)c(N)c(Cl)c1",         # clenbuterol
        "NC(CC(C)(C)C)(C)C",                        # 2,4,4-trimethyl-2-Pentanamine
        "CC(C)(C)NC[C@@H](O)c1cc(Cl)c(N)c(Cl)c1",    # (S)-clenbuterol
        "NC1CC1",                                  # cyclopropylamine
        "C[C@H](N)c1ccccc1",                        # (1S)-1-phenylethanamine
        "Nc1ccc(cc1)\\N=N\\c1ccc(N)cc1",             # 4,4'-diaminoazobenzene
        "CO\\C=C\\C(=O)C1=NC2=C3C(C=NC3=C(N)C=C2NC(C)=O)=C1",  # lymphostin
        "NCCCCCCCCCCC[C@H](CC)C",                   # Medelamine B
        "Nc1cnnc(O)c1Cl",                          # chloridazone-desphenyl
        "Cl.Cl.NCCCCN",                            # 1,4-Diaminobutane dihydrochloride
        "NCC1(CCCC1)c1ccccc1",                      # 1-(1-phenylcyclopentyl)methylamine
        "C[C@H](C(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)O)NC(=O)[C@@H](CC3=CC=CC=C3)N",  # sevadicin
        "COC(=O)[C@H]1[C@@H]([C@H]2CN3C(=O)C=CC=C3[C@@H]1N2CCC4=CC=CC=C4)CO",  # LSM-13982
        "NCCC=C",                                  # 3-buten-1-amine
        "CC(=O)[C@H](CCCCN)NS(=O)(=O)c1ccc(C)cc1",  # p-Ts-L-Lys-Me
        "NCc1ccccc1",                              # benzylamine
        # Additional test SMILES can be added here.
    ]
    
    for smi in test_smiles:
        result, reason = is_primary_amine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)