"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: Carboxamidine-containing compounds.
Definition: Compounds having the structure RC(=NR)NR2,
i.e. the –C(=NH)NH2 group or a substituted analogue is present.
This version refines matching by:
  • Using two SMARTS patterns (neutral and cationic)
  • Checking that the double‐bonded N is not directly bound to an –OH moiety
    (which would indicate an amidoxime rather than a carboxamidine)
  • Looking at the extra substituent on the carboxamidine carbon to rule out
    cases where it is bonded to an unexpected heteroatom (as in azo or other systems)
  • Relaxing the peptide filter: only very heavy molecules with many amide bonds are rejected.
  
Examples provided by the user include formamidine, acetamidine, benzamidine, guanidine etc.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group (RC(=NR)NR2).

    The algorithm:
      1. Parses the input SMILES.
      2. Adds hydrogen explicitly (to allow neighbor analysis).
      3. Uses two SMARTS patterns to capture both neutral and cationic representations.
         The SMARTS pattern "[CX3](=[NX2])[NX3]" matches a trigonal carbon double bonded to a nitrogen and
         singly bonded to another nitrogen.
      4. For each match, check:
         - That the double-bonded nitrogen is not bound (by a single bond) to any oxygen that carries at least one hydrogen.
           (This helps to filter out amidoxime-like groups.)
         - That the central carbon has one extra neighbor (its R group) and if present, that the extra substituent
           is acceptable (either carbon or hydrogen). If the extra neighbor is not present (as in formamidine) that is acceptable.
      5. Counts amide bonds (using a simple "C(=O)N" SMARTS) so that if a very heavy molecule (mol weight >600) shows
         many (>=4) amide bonds it is flagged as peptide-like (and rejected).
    Returns:
      bool: True if the molecule is classified as containing a carboxamidine group.
      str: Explanation for the classification.
    """
    # Parse SMILES and add explicit hydrogens for correct neighbor analysis.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    mol = Chem.AddHs(mol)
    
    # Define two SMARTS patterns to capture the carboxamidine motif.
    # These patterns look for a carbon with a double bond to a nitrogen and a single bond to another nitrogen.
    pattern_neutral = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    pattern_cationic = Chem.MolFromSmarts("[CX3](=[NX2+])[NX3]")
    if pattern_neutral is None or pattern_cationic is None:
        return False, "Error defining SMARTS patterns."
    
    # Get all matches from both patterns.
    matches1 = mol.GetSubstructMatches(pattern_neutral)
    matches2 = mol.GetSubstructMatches(pattern_cationic)
    # Use a set of tuples to avoid duplicates.
    all_matches = set(matches1) | set(matches2)
    if not all_matches:
        return False, "Carboxamidine moiety (-C(=NH)NH2 or substituted equivalent) not found."
    
    valid_matches = []
    for match in all_matches:
        # Expect match order: (idx_c, idx_dblN, idx_singleN)
        idx_c, idx_ndbl, idx_nsingle = match
        atom_c = mol.GetAtomWithIdx(idx_c)
        atom_ndbl = mol.GetAtomWithIdx(idx_ndbl)
        atom_nsingle = mol.GetAtomWithIdx(idx_nsingle)
        
        # ------------------------------
        # (A) Filter based on oxygen on the double-bonded N:
        # If the double bonded N has any neighbor oxygen where the bond type is SINGLE and that oxygen carries at least one hydrogen,
        # then skip this match as it likely is an amidoxime.
        skip_match = False
        for nbr in atom_ndbl.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in (idx_c, idx_nsingle):
                # Get the bond between atom_ndbl and the oxygen
                bond = mol.GetBondBetweenAtoms(atom_ndbl.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    # Check if oxygen has any explicit hydrogen attached.
                    if nbr.GetTotalNumHs() > 0:
                        skip_match = True
                        break
        if skip_match:
            continue
        
        # ------------------------------
        # (B) Check the extra substituent on the carboxamidine carbon.
        # The carbon is already bonded to the double-bonded N and the single-bonded N.
        # It may have a third neighbor (its R group). For formamidine there is no extra substituent (implicit hydrogen is fine).
        # If an extra neighbor exists, we require it to be a carbon (or hydrogen) rather than an unexpected heteroatom.
        extra_neighbours = [nbr for nbr in atom_c.GetNeighbors() if nbr.GetIdx() not in (idx_ndbl, idx_nsingle)]
        if extra_neighbours:
            for nbr in extra_neighbours:
                if nbr.GetAtomicNum() not in (1, 6):  # allow hydrogen or carbon.
                    # In some systems the substituent may be part of an aromatic ring but then it should be a carbon.
                    skip_match = True
                    break
        if skip_match:
            continue
        
        valid_matches.append(match)
    
    if not valid_matches:
        return False, "Carboxamidine moiety not found after refined filtering."
    
    # Count amide bonds in the molecule (using a simple SMARTS for C(=O)N).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # If the molecule is very heavy and has many amide bonds, it is likely peptide‐like.
    if num_amide >= 4 and mol_wt > 600:
        return False, "Molecule appears to be peptide-like (many amide bonds in a heavy molecule), likely a false positive."
    
    return True, f"Found carboxamidine group in {len(valid_matches)} location(s)."

# For testing/debugging only (run when executed as a script).
if __name__ == "__main__":
    test_smiles = [
        "[H]C(N)=N",  # formamidine: expected true positive.
        "CN(Cc1ccc(Cl)nc1)C(\\C)=N/C#N",  # (Z)-acetamiprid: true positive.
        "CC(=N)NCC1=CC(CN)=CC=C1",       # N-[3-(aminomethyl)benzyl]acetamidine: true positive.
        "CN1CCCN=C1\\C=C\\c1cccs1",       # pyrantel: true positive.
        "CC1=N[C@@H]([C@@H](O)CN1)C(O)=O",# 5-hydroxyectoine: true positive.
        "CN(Cc1ccc(Cl)nc1)C(C)=NC#N",      # acetamiprid: true positive.
        "[H][C@@]1(COc2ccc(cc2)-c2ccc(cc2)C(N)=N)C[C@@H](CC(O)=O)C(=O)N1",  # fradafiban: true positive.
        "P(=O)(N=C(N(CC)CC)C)(OC)F",      # A-232 nerve agent: true positive.
        "NC(=N)c1ccc2cccc(Nc3ncccn3)c2c1",# 8-(pyrimidin-2-ylamino)naphthalene-2-carboximidamide: true positive.
        "CC(N)=N",                      # acetamidine: true positive.
        "O[C@H]1[C@@H](O)C(NC(=N)CNC=O)O[C@@H]1COP(O)(O)=O",  # 2-formamido-N(1)-(5-phospho-D-ribosyl)acetamidine: true positive.
        "N[C@@H](CCCCC(N)=N)C(O)=O",      # L-indospicine: true positive.
        "NC(=N)N1CC(O)c2ccccc2C1",        # 4-hydroxydebrisoquin: true positive.
        "Cc1cc(c(O)c(C)c1CC1=NCCN1)C(C)(C)C",  # oxymetazoline: true positive.
        "CC1([C@@H](N2[C@@H](S1)[C@@H](C2=O)N=CN3CCCCCC3)C(=O)OCOC(=O)C(C)(C)C)C",  # 2,2-dimethylpropanoyloxymethyl ...: true positive.
        "CN(C=Nc1ccc(C)cc1C)C=Nc1ccc(C)cc1C",# amitraz: true positive.
        "CC(Oc1c(Cl)cccc1Cl)C1=NCCN1",      # lofexidine: true positive.
        "CN=CNc1ccc(C)cc1C",              # N'-(2,4-dimethylphenyl)-N-methylformamidine: true positive.
        "NC(=N)N1CCc2ccccc2C1",           # debrisoquin: true positive.
        "P(=O)(N=C(N(CC)CC)C)(C)F",        # A-230 nerve agent: true positive.
        "CC1=N[C@@H](CCN1)C(O)=O",         # ectoine: true positive.
        "CCCCCNC(=N)N\\N=C\\c1c[nH]c2ccc(OC)cc12",  # tegaserod: true positive.
        "NC(=N)c1ccc(NN=Nc2ccc(cc2)C(N)=N)cc1",  # diminazene: true positive.
        "NC(=N)c1c[nH]c2nc(N)[nH]c(=O)c12",      # 7-formamidino-7-deazaguanine: true positive.
        "C(CCCOC1=CC=C(C=C1)C(N)=N)COC2=CC=C(C=C2)C(N)=N",  # pentamidine: true positive.
        "C[C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]2O)[C@@H](N)C[C@@H]1NC(=N)C(O)=O",  # kasugamycin: true positive.
        "CON=CNC(=O)c1ccc(cc1C)C1=NO[C@@](C1)(c1cc(Cl)cc(Cl)c1)C(F)(F)F",  # (S)-fluxametamide: true positive expected but may be missed if too many oxygen neighbors occur.
        "CN(C)C=Nc1cccc(O)c1",  # N'-(3-hydroxyphenyl)-N,N-dimethylformamidine: true positive.
        "CC\\N=C(/N)c1ccc(cc1)-c1ccc(o1)-c1ccc(cc1)C(\\N)=N\\CC",  # 2,5-bis{[4-(N-ethylamidino)]phenyl}furan: true positive.
        # False positives examples (these should NOT be classified as carboxamidine):
        "CC(C)(C(=N)N)N=NC(C)(C)C(=N)N",  # wrongly classified: extra unsaturated groups etc.
        "[H+].[H+].[O-]C(=O)\\C=C/C([O-])=O.CCCCCNC(=N)N\\N=C\\c1c[nH]c2ccc(OC)cc12",  # tegaserod maleate: should be rejected.
        "C1(=NC=C(C=C1C(=O)O)COC)C2=NC(C(N2)(C)C(C)C)=O",  # another false positive.
        "C1C(N=C(C(=O)O1)NNC2=CC=CC=C2Cl)(CO)CO",  # another false positive.
        "CC1=CC=C(C=C1)C2=NC(C3=C(N2)N(C(=O)[C@H](NCC(O)=O)C1CCCCC1)C4=CC=C(C=C4)F)(C(F)(F)F)C(F)(F)F",  # another false positive.
        "CC1=CC=C(C=C1)N(CC2=NCCN2)C3=CC=C(C=C3)O",  # another false positive.
        "C=1(N(C(N=C(N1)N)(C)C)C=2C=CC(=CC2)Cl)N.C3=CC=CC4=C3C=C(C(=C4CC=5C6=C(C=C(C5O)C(=O)O)C=CC=C6)O)C(=O)O.C=7(N(C(N=C(N7)N)(C)C)C=8C=CC(=CC8)Cl)N",  # Cycloguanil pamoate: false positive.
        "S(=O)(=O)(O[C@@H]1C(O)(O)[C@]23NC(N)=N[C@H]2[C@H](CO)N=C(N3C1)N",  # Decarbamoylgonyautoxin-3: false positive.
        "OC(=O)CNC=N",  # N-formimidoylglycine: false positive.
        "S(/C(=N\\C#N)/N1CCC(C(=O)N)CC1)C",  # Methyl 4-(aminocarbonyl)-n-cyanotetrahydropyridine-1(2H)-carbimidothioate: false positive.
        "NC(=N)NCCCCCC(O)=O",  # 6-guanidinohexanoic acid: false positive.
        "O.Cl.Cl.C\\C(\\C=N\\NC(N)=N)=N/NC(N)=N",  # mitoguazone hydrochloride hydrate: false positive.
        "O(C1C(O)C(O)C(OC1)(O)CNC(CCCN=C(N)N)C(O)=O)C2OC(C(O)C(O)C2O)CO",  # N2-Maltulosylarginine: false positive.
        "CN1CCN(CC1)C=NC2=C(C3=C(S2)CCCCC3)C#N",  # another false positive.
        "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1N=CNCC2=O",  # dehydrocoformycin: false positive.
        "OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H]",  # L-arginine-d7: false positive.
        "CN1COCN(Cc2cnc(Cl)s2)C1=N[N+]([O-])=O",  # thiamethoxam: false positive.
        "CC1=CC=C(C=C1)NC(=O)CSC2=NCCN2",  # another false positive.
        "C1CCC2(C1)C3=C(CCCC3=O)NC(=N2)NC4=NC5=CC=CC=C5O4",  # another false positive.
        "CC1=CC(=CC=C1)NC(=S)N=C(N)NC2=NC3=C(C=CC(=C3)C)C(=N2)C",  # another false positive.
        "COC1=CC(=CC(=C1O)OC)C=C2C(=O)N=C(S2)NC3=CC=C(C=C3)O",  # another false positive.
        "C1=CC=C(C=C1)CNC2=NC(=O)C(=CC3=CC=CC=C3)S2",  # another false positive.
        "COC1=C(C=CC(=C1)C=C2C(=O)N=C(S2)NC3=CC=CC=C3Cl)O",  # another false positive.
        "C1CNC2=NCCCN2C1"  # another false positive.
    ]
    
    for smi in test_smiles:
        res, reason = is_carboxamidine(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")