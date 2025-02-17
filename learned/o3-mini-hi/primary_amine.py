"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary amine
A primary amine is defined as a compound derived from ammonia by the substitution of exactly one hydrogen 
atom by a hydrocarbyl group. In practice, we look for an –NH2 group (a trivalent nitrogen with exactly two hydrogens)
attached to exactly one heavy (non‐hydrogen) atom, and we also exclude cases where that heavy neighbor is part
of a carbonyl (C=O), which would indicate an amide.
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule contains at least one primary amine group based on its SMILES string.
    A primary amine group is defined as an NH2 group attached to exactly one heavy (non-hydrogen) atom;
    that heavy neighbor must not be a carbon that is double-bonded to oxygen (i.e. not part of an amide).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one primary amine group is detected, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens so that we have both explicit and implicit counts available.
    mol = Chem.AddHs(mol)
    
    # Loop over each nitrogen atom.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Only interested in nitrogen atoms.
        
        # Instead of GetTotalNumHs, explicitly sum explicit and implicit hydrogens.
        num_H = atom.GetNumExplicitHs() + atom.GetNumImplicitHs()
        # For a primary amine (R–NH2), expect exactly two hydrogens on the nitrogen.
        if num_H != 2:
            continue
        
        # Determine non-hydrogen (heavy) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            continue  # Must have exactly one heavy atom neighbor.
        
        # Check that the heavy neighbor is not part of a carbonyl group (i.e. not adjacent to C=O).
        heavy_nbr = heavy_neighbors[0]
        carbonyl_found = False
        if heavy_nbr.GetAtomicNum() == 6:  # Only carbon can be the carbonyl carbon.
            for bond in heavy_nbr.GetBonds():
                # Look for a double bond to oxygen.
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(heavy_nbr)
                    if other.GetAtomicNum() == 8:
                        carbonyl_found = True
                        break
        if carbonyl_found:
            continue
        
        # If we reach here for any nitrogen atom, we have found a valid primary amine group.
        return True, "Contains a primary amine group (R–NH2): an NH2 group attached to exactly one heavy atom without adjacent carbonyl."
    
    return False, "No primary amine group (R–NH2) found"

# Example usage (for testing a few cases)
if __name__ == "__main__":
    test_smiles = [
        "Nc1ccc(cc1)\\N=N\\c1ccccc1",       # 4-(phenylazo)aniline
        "C[C@@H](N)Cc1ccccc1",               # (R)-amphetamine
        "C1=CC=CC2=C1C(C3=C(C2=O)C(=C(C(=C3N)O)S(O)(=O)=O)O)=O",  # nuclear fast red free acid
        "CC(C)NCC(C)(C)N",                  # N(1)-isopropyl-2-methylpropan-1,2-diamine
        "CCCCCCCCCCCCCCCCCCN",              # octadecan-1-amine
        "[H]C(=O)CCCNCCCN",                 # N-(3-aminopropyl)-4-aminobutanal
        "CC(C)C(N)C(C)C",                  # 2,4-dimethylpentan-3-amine
        "CCCCN",                           # butan-1-amine
        "C1=CC(=C(C=C1O)N)CCN",             # 3-amino-4-(2-aminoethyl)phenol
        "Nc1ccc2ccccc2c1",                 # 2-naphthylamine
        "CC(C)(N)Cc1ccccc1",                # phentermine
        "CN",                              # methylamine
        "COC1=C(O)C=C(CCN)C=C1",            # 4-methoxytyramine
        "CC(C)(C)NCC(O)c1cc(Cl)c(N)c(Cl)c1", # clenbuterol
        "NC(CC(C)(C)C)(C)C",               # 2,4,4-trimethyl-2-Pentanamine
        "CC(C)(C)NC[C@@H](O)c1cc(Cl)c(N)c(Cl)c1",  # (S)-clenbuterol
        "NC1CC1",                         # cyclopropylamine
        "C[C@H](N)c1ccccc1",               # (1S)-1-phenylethanamine
        "Nc1ccc(cc1)\\N=N\\c1ccc(N)cc1",    # 4,4'-diaminoazobenzene
        "CO\\C=C\\C(=O)C1=NC2=C3C(C=NC3=C(N)C=C2NC(C)=O)=C1",  # lymphostin
        "NCCCCCCCCCCC[C@H](CC)C",          # Medelamine B
        "Nc1cnnc(O)c1Cl",                 # chloridazone-desphenyl
        "Cl.Cl.NCCCCN",                   # 1,4-Diaminobutane dihydrochloride
        "NCC1(CCCC1)c1ccccc1",             # 1-(1-phenylcyclopentyl)methylamine
        "C[C@H](C(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)O)NC(=O)[C@@H](CC3=CC=CC=C3)N",  # sevadicin
        "COC(=O)[C@H]1[C@@H]([C@H]2CN3C(=O)C=CC=C3[C@@H]1N2CCC4=CC=CC=C4)CO",  # LSM-13982
        "NCCC=C",                         # 3-buten-1-amine
        "CC(=O)[C@H](CCCCN)NS(=O)(=O)c1ccc(C)cc1",  # p-Ts-L-Lys-Me
        "NCc1ccccc1",                     # benzylamine
        # ... further test SMILES can be added here.
    ]
    
    for smi in test_smiles:
        result, reason = is_primary_amine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)