"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary amine
A primary amine is formally derived from ammonia (NH3) by replacing exactly one hydrogen with a hydrocarbyl group.
Thus the defining functional group is –NH2 connected to exactly one non‐hydrogen atom (usually carbon) and not adjacent to a carbonyl (i.e. not an amide).
In addition, if a molecule is large and contains several amide bonds (as in a peptide), then the incidental free NH2 is not considered to define
the chemical class “primary amine”.
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is (primarily) a primary amine based on its SMILES string.
    
    The algorithm works by:
      1. Parsing the SMILES string and adding explicit hydrogens.
      2. Searching for at least one nitrogen atom that has exactly 2 hydrogens (when counted with GetTotalNumHs)
         and exactly one heavy-atom (non-H) neighbor.
      3. Discarding candidate N if its single heavy neighbor is a carbon that is double-bonded to an oxygen (i.e. part of an amide group).
      4. As a heuristic, if the molecule is large (many heavy atoms) and it has more than one amide bond (substructure C(=O)N),
         we assume it is a multifunctional (e.g. peptide) compound and not classified as simply a primary amine.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is considered a primary amine, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are correct.
    mol = Chem.AddHs(mol)
    
    # Count amide bonds (pattern: C(=O)N) in the molecule.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    num_amide = len(mol.GetSubstructMatches(amide_smarts))
    
    primary_count = 0  # count of nitrogen atoms that qualify as primary amine groups
    
    for atom in mol.GetAtoms():
        # Focus only on nitrogen atoms.
        if atom.GetAtomicNum() != 7:
            continue
        # Get total number of hydrogens attached (explicit+implicit now explicit).
        num_H = atom.GetTotalNumHs()
        # A primary amine (R–NH2) must have exactly two hydrogens.
        if num_H != 2:
            continue
        
        # Get heavy-atom neighbors (non-hydrogen connected atoms)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # For a primary amine group the nitrogen should have exactly one heavy-atom neighbor.
        if len(heavy_neighbors) != 1:
            continue
        
        neighbor = heavy_neighbors[0]
        # Exclude cases where the nitrogen is attached to a carbonyl group (i.e. an amide).
        is_adjacent_to_carbonyl = False
        for bond in neighbor.GetBonds():
            # Check if the bond is a double bond.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(neighbor)
                # If the neighbor of the heavy atom is oxygen then assume a carbonyl.
                if other.GetAtomicNum() == 8:
                    is_adjacent_to_carbonyl = True
                    break
        if is_adjacent_to_carbonyl:
            continue
        
        # If we reach here, we have found a candidate primary amine group.
        primary_count += 1

    # In case no qualifying primary amine group was found:
    if primary_count == 0:
        return False, "No primary amine group (R–NH2) found"
    
    # Heuristic to avoid false positives in large multifunctional molecules (e.g. peptides)
    # If the molecule is large (many heavy atoms) and has more than one amide bond, it might be a peptide
    # in which the free amine is incidental.
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_amide > 1 and num_heavy_atoms > 50:
        return False, ("Molecule appears to be large and contains multiple amide bonds, " 
                       "suggesting it is a peptide or multifunctional compound rather than a primary amine")
    
    return True, "Contains a primary amine group (R–NH2) as defined by exactly one non‐hydrogen substituent and two hydrogens"

    
# Example usage for testing on a few SMILES (feel free to add more cases):
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
        "COC(=O)[C@H]1[C@@H]([C@H]2CN3C(=O)C=CC=C3[C@@H]1N2CCC4=CC=CC=C4)CO",  # LSM-13982 (note: similar to LSM-12804)
        "NCCC=C",                                  # 3-buten-1-amine
        "CC(=O)[C@H](CCCCN)NS(=O)(=O)c1ccc(C)cc1",  # p-Ts-L-Lys-Me
        "NCc1ccccc1",                              # benzylamine
        # ... additional test SMILES can be added here.
    ]
    
    for smi in test_smiles:
        result, reason = is_primary_amine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)