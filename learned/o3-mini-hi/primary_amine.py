"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary amine
A primary amine is a compound formally derived from ammonia (NH3) by replacing exactly one hydrogen 
with a hydrocarbyl group. Thus the key functional group is –NH2 (with exactly one non‐hydrogen substituent)
that is not adjacent to a carbonyl group (i.e. not part of an amide).
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule contains at least one primary amine group based on its SMILES string.
    
    A primary amine group is defined as an -NH2 group attached to exactly one heavy (non-hydrogen) atom,
    with that neighbor not being part of a carbonyl (C=O) bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one primary amine group is detected, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to reliably count them.
    mol = Chem.AddHs(mol)
    
    found_primary = False  # Flag to indicate a primary amine group was found
    
    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Skip non-nitrogen atoms
        
        # Get the total number of hydrogen atoms attached to this nitrogen.
        num_h = atom.GetTotalNumHs()
        # For a primary amine (R–NH2), expect exactly two hydrogens.
        if num_h != 2:
            continue
        
        # Count heavy (non hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            continue  # Must have exactly one non-hydrogen substituent
        
        # Check that the heavy substituent is not part of a carbonyl group.
        neighbor = heavy_neighbors[0]
        is_carbonyl = False
        # If the neighbor is carbon, check if it is double bonded to oxygen.
        if neighbor.GetAtomicNum() == 6:
            for bond in neighbor.GetBonds():
                # Ensure we are looking at bonds where neighbor is involved.
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(neighbor)
                    if other.GetAtomicNum() == 8:
                        is_carbonyl = True
                        break
        if is_carbonyl:
            continue  # Exclude amide-like connectivity
        
        # If we reached here then we found a candidate primary amine (-NH2) group.
        found_primary = True
        break  # No need to check further
    
    if found_primary:
        return True, "Contains at least one primary amine group (R–NH2), defined as an -NH2 group attached to exactly one heavy atom and with no adjacent carbonyl."
    else:
        return False, "No primary amine group (R–NH2) found"
    
# Example usage (for testing a few cases)
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
        # Further test SMILES can be added here
    ]
    
    for smi in test_smiles:
        result, reason = is_primary_amine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)