"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary amine
A primary amine is defined as a compound formally derived from ammonia by replacing one hydrogen atom 
with a hydrocarbyl group. In practice, we look for an –NH2 group (a trivalent nitrogen bearing exactly two hydrogens)
that is attached to exactly one heavy (non‐hydrogen) atom, excluding cases where that heavy neighbor is part 
of a carbonyl group (C=O), which would indicate an amide.
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule contains at least one primary amine group (R-NH2) based on its SMILES string.
    The group is identified as a nitrogen atom with exactly two attached hydrogens that is not adjacent to a 
    carbonyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one primary amine group is detected, False otherwise.
        str : Explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a primary amine:
    # [NX3;H2] : A trivalent nitrogen with exactly 2 hydrogens.
    # !$(NC=O)  : Excludes cases where the nitrogen is attached to a carbonyl group.
    primary_amine_smarts = "[NX3;H2;!$(NC=O)]"
    patt = Chem.MolFromSmarts(primary_amine_smarts)
    if patt is None:
        return False, "Error in constructing SMARTS pattern"
    
    # Find substructure matches in the molecule.
    matches = mol.GetSubstructMatches(patt)
    if matches:
        # If at least one match found, we conclude that the molecule contains a primary amine.
        return True, "Contains a primary amine group (R–NH2): a nitrogen with exactly two hydrogens not adjacent to a carbonyl."
    else:
        return False, "No primary amine group (R–NH2) found"

# Example usage for testing various cases:
if __name__ == "__main__":
    test_cases = [
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
        "[O-]/[N+](/N(CCCN)CCCN)=N\\O",     # N-[Bis(3-aminopropyl)amino]-N-hydroxynitrous amide
        "CC1=C(ON=C1C2CCN(CC2)CCC3=CC=CC=C3)C4=CC=C(C=C4)Cl",  # Complex example
        "C1=CC(=CC=C1)CC2(C3=C(C(N4C2(CCC4=O)O)C(NCCC(=O)O)=O)C=C(C=C3O)O)N",  # perquinoline B
        "C1CCC(C(C1)N)N",                 # cyclohexane-1,2-diamine
        "NC1CCC(CCC)CC1",                 # 4-Propylcyclohexylamine
        "CC(C)(C)NC[C@H](O)c1cc(Cl)c(N)c(Cl)c1",  # (R)-clenbuterol
        "Cc1cc(ccc1N=Nc1ccc2c(cc(c(N)c2c1O)S(O)(=O)=O)S(O)(=O)=O)-c1ccc(N=Nc2ccc3c(cc(c(N)c3c2O)S(O)(=O)=O)S(O)(=O)=O)c(C)c1",  # Evans blue free acid
        "N(C(CCCN)([H])[H])([H])[H]",    # Putrescine_d4
        "C1=CC(=C(C(=C1)O)O)CCN",          # 3-(2-aminoethyl)benzene-1,2-diol
        "NC[C@H](CC)C",                  # 2-Methylbutylamine
        "Nc1ccccc1",                     # aniline
        "[O-]/[N+](/N(CCCN)CCC)=N\\O",     # PAPA NONOate
        "C1(=CC(=CC=C1O)CCN)OC",           # 3-methoxytyramine
        "NCCCCCCCCCCCCCC",                # Tetradecylamine
        "CN[C@@H]1[C@H](O)[C@H](NC)[C@H]2O[C@]3(O)[C@@H](O[C@H](C)C=C3OC)O[C@@H]2[C@H]1O",  # Spenolimycin
        "CC(C)=CCCN",                    # 4-methyl-3-penten-1-amine
        "Cc1cccc(C)c1N",                 # 2,6-dimethylaniline
        "Nc1cc2ccccc2c2ccccc12"           # phenanthren-9-amine
    ]
    
    for smi in test_cases:
        result, reason = is_primary_amine(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)