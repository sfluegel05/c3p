"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: Gas Molecular Entity
Definition: Any main group molecular entity that is gaseous at STP (0°C and 100 kPa).

This heuristic classifier works roughly as follows:
  1. Parse the SMILES and add explicit hydrogens.
  2. If the molecule is monoatomic, only noble gases are accepted.
  3. For polyatomic species the number of heavy atoms (i.e. non‐H atoms) must be small (≤4)
     and all heavy atoms must be “safe” (only C, N, O, F, Cl, Br, I).
  4. However, two well–known exceptions are allowed even if formal charges are present:
       – carbon monoxide ([C-]#[O+]) and
       – ozone (all O’s with the canonical O–O+–O pattern).
  5. Many “false positives” (such as methanediimine, chloromethane, methyloxidanyl, etc.) are
     rejected by further requiring that:
       • For diatomic molecules (n_heavy ==2) the only allowed heteronuclear combinations are those
         that occur in CO or for a halogen diatomic (which must be homonuclear) or, for C–Cl,
         only when the Cl is in a vinyl group (e.g. chloroethene).
       • Molecules composed solely of C and N (with 3 or 4 heavy atoms) tend to be rejected.
       • Molecules containing oxygen are accepted only if they occur in CO, CO2, oxirane or ozone.
       • Certain polar functional groups (–OH, –CHO, –COOH) are rejected.
       • A few forbidden motifs (such as a C≡N triple bond or isocyanate) are looked for.
  6. If all tests are passed, the molecule is classified as a gas molecular entity.
    
Note: This is a heuristic classifier using structural surrogates for volatility.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule (as a SMILES string) is a main group gas molecular entity at STP.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple of (True if accepted, False otherwise, explanation string)
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that all atoms are counted.
    mol = Chem.AddHs(mol)
    
    # Require that the molecule is a single connected entity.
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Molecule is disconnected; should be a single connected entity"
    
    total_atoms = mol.GetNumAtoms()  # includes explicit Hs
    # Get heavy atoms, i.e. atoms with atomic number != 1
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    n_heavy = len(heavy_atoms)
    heavy_symbols = [atom.GetSymbol() for atom in heavy_atoms]
    
    # Special case: monoatomic species (note: without added H, monoatomic SMILES have 1 atom).
    if total_atoms == 1:
        # Only accept if the atom is a noble gas.
        noble_gases = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}
        if heavy_symbols[0] in noble_gases:
            mw = Descriptors.ExactMolWt(mol)
            return True, f"Monoatomic noble gas ({heavy_symbols[0]}, MW {mw:.1f} Da) is classified as a gas"
        else:
            return False, f"Monoatomic {heavy_symbols[0]} is not classified as a gas (only noble gases accepted)"
    
    # Identify two known exceptions that are allowed even if formal charges are present.
    # Exception: carbon monoxide.
    is_CO = (n_heavy == 2 and sorted(heavy_symbols) == ["C", "O"])
    # Exception: ozone – three oxygens.
    is_O3 = (n_heavy == 3 and all(sym == "O" for sym in heavy_symbols))
    
    # Check that each heavy atom has formal charge zero (unless the molecule is CO or O3).
    if not (is_CO or is_O3):
        for atom in heavy_atoms:
            if atom.GetFormalCharge() != 0:
                return False, f"Atom {atom.GetSymbol()} carries nonzero formal charge; molecule rejected"
    
    # Enforce that typical gas molecular entities have at most 4 heavy atoms.
    if n_heavy > 4:
        return False, f"Too many heavy atoms ({n_heavy}); typical gas molecular entities have ≤4 heavy atoms"
    
    hs_set = set(heavy_symbols)
    
    # ---- Additional ad hoc rejection rules ----
    # For diatomic (n_heavy==2) molecules:
    if n_heavy == 2:
        # If both heavy atoms are halogens, allow only if they are identical (e.g. FF, ClCl).
        halogens = {"F", "Cl", "Br", "I"}
        if hs_set.issubset(halogens):
            if heavy_symbols[0] != heavy_symbols[1]:
                return False, f"Heteronuclear diatomic halogen ({heavy_symbols[0]}, {heavy_symbols[1]}) not accepted as a gas"
        # For a pair of C and O heavy atoms, only allow carbon monoxide.
        if hs_set == {"C", "O"} and not is_CO:
            return False, "Molecule with C and O heavy atoms (likely formaldehyde) is rejected"
        # For a C–Cl combination, allow only if the Cl is bonded to an sp2 carbon (as in vinyl chloride).
        if hs_set == {"C", "Cl"}:
            pat_vinyl = Chem.MolFromSmarts("[C]=[C]-Cl")
            if not (pat_vinyl and mol.HasSubstructMatch(pat_vinyl)):
                return False, "Molecule with C and Cl heavy atoms not matching a vinyl chloride motif; likely not a gas"
    
    # For molecules composed solely of C and N: if there are 3 or 4 heavy atoms, reject (e.g. methanediimine, trimethylamine)
    if hs_set == {"C", "N"} and n_heavy in {3, 4}:
        return False, f"Molecule with heavy atoms {sorted(heavy_symbols)} (e.g. methanediimine or trimethylamine) rejected"
    
    # Reject molecules containing polar functional groups that would tend to raise the boiling point.
    pat_oh = Chem.MolFromSmarts("[OX2H]")  # hydroxyl group
    if mol.HasSubstructMatch(pat_oh):
        return False, "Molecule contains an -OH group; likely not a gas"
    pat_aldehyde = Chem.MolFromSmarts("[CX3H1](=O)")
    if mol.HasSubstructMatch(pat_aldehyde):
        return False, "Molecule contains an aldehyde group; likely not a gas"
    pat_carboxylic = Chem.MolFromSmarts("C(=O)[OH]")
    if mol.HasSubstructMatch(pat_carboxylic):
        return False, "Molecule contains a carboxylic acid group; likely not a gas"
    
    # Look for forbidden bonding motifs.
    # Reject if any triple bond occurs between a carbon and a nitrogen (nitrile motif).
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            at1 = bond.GetBeginAtom().GetSymbol()
            at2 = bond.GetEndAtom().GetSymbol()
            if {"C", "N"} == {at1, at2}:
                return False, "Molecule has a C≡N triple bond (nitrile motif); likely not a gas"
    # Reject isocyanate motifs.
    pat_iso = Chem.MolFromSmarts("[N]=C=O")
    if pat_iso and mol.HasSubstructMatch(pat_iso):
        return False, "Molecule contains isocyanate motif [N]=C=O; likely not a gas"
    
    # If the molecule contains oxygen, only allow if it fits an allowed pattern:
    # allowed patterns include carbon dioxide, oxirane, CO (already handled) or ozone.
    if "O" in heavy_symbols:
        pat_co2 = Chem.MolFromSmarts("O=C=O")
        pat_oxirane = Chem.MolFromSmarts("C1CO1")
        if not (mol.HasSubstructMatch(pat_co2) or mol.HasSubstructMatch(pat_oxirane) or is_CO or is_O3):
            return False, "Molecule contains oxygen outside allowed patterns; likely not a gas"
    
    # If we pass all tests, accept the molecule.
    mw = Descriptors.ExactMolWt(mol)
    halogen_nums = {9, 17, 35, 53}
    n_halogen = sum(1 for a in heavy_atoms if a.GetAtomicNum() in halogen_nums)
    halogen_fraction = n_halogen / n_heavy if n_heavy else 0
    return True, (f"Molecule (MW {mw:.1f} Da, {n_heavy} heavy atoms, halogen fraction {halogen_fraction:.2f}) "
                  "is classified as a gas molecular entity")

# Example usage:
if __name__ == '__main__':
    test_cases = [
        # True positives
        ("[H]\\C(C)=C(\\[H])C", "trans-but-2-ene"),
        ("FC=C", "fluoroethene"),
        ("[H]N([H])[H]", "ammonia"),
        ("Cl[H]", "hydrogen chloride"),
        ("[H]C([H])([H])[H]", "methane"),
        ("[3He]", "helium-3 atom"),
        ("C=C", "ethene"),
        ("[1H][1H]", "diprotium"),
        ("[220Rn]", "radon-220 atom"),
        ("[Rn]", "radon(0)"),
        ("[H][H]", "dihydrogen"),
        ("[He]", "helium(0)"),
        ("FF", "difluorine"),
        ("[6He]", "helium-6 atom"),
        ("CC", "ethane"),
        ("ClC=C", "chloroethene"),
        ("[219Rn]", "radon-219 atom"),
        ("ClCl", "dichlorine"),
        ("[Xe]", "xenon atom"),
        ("[H]\\C(C)=C(/[H])C", "cis-but-2-ene"),
        ("I[H]", "hydrogen iodide"),
        ("CCC=C", "but-1-ene"),
        ("CC#C", "propyne"),
        ("CC(C)C", "isobutane"),
        ("[H]C(C)=C([H])C", "but-2-ene"),
        ("[O][O]", "triplet dioxygen"),
        ("[4He]", "helium-4 atom"),
        ("[Kr]", "krypton atom"),
        ("C1CO1", "oxirane"),
        ("[Ar]", "argon atom"),
        ("O=C=O", "carbon dioxide"),
        ("[222Rn]", "radon-222 atom"),
        ("O=[13C]=O", "((13)C)carbon dioxide"),
        ("[3H][3H]", "ditritium"),
        ("[Ne]", "neon(0)"),
        ("CCCC", "butane"),
        ("CCC", "propane"),
        ("CC(C)=C", "2-methylprop-1-ene"),
        # False negatives (should be accepted)
        ("[C-]#[O+]", "carbon monoxide"),
        # Borderline cases (many of these will be rejected):
        ("N=C=N", "methanediimine"),
        ("[H]C([H])([H])Cl", "chloromethane"),
        ("C[O]", "methyloxidanyl"),
        ("[H][C]=O", "oxomethyl"),
        ("NC(O)=O", "carbamic acid"),
        ("[H][C@@]1(C)CO1", "(R)-1,2-epoxypropane"),
        ("C1CC1O", "cyclopropanol"),
        ("[Xe+]", "xenon(.1+)"),
        ("[C][H]", "hydridocarbon(.)"),
        ("[1H][1H+]", "diprotium(.1+)"),
        ("O(N)C", "Methoxyamine"),
        ("C1OC=C1", "oxetene"),
        ("[He++]", "helium(2+)"),
        ("C1N=N1", "3H-diazirine"),
        ("O(C(=O)C([2H])([2H])[2H])[2H]", "acetic acid-d4"),
        ("C[C]=O", "acetyl"),
        ("[H][N]O[H]", "hydroxyazanyl"),
        ("Cl\\C=C\\Cl", "trans-1,2-dichloroethene"),
        ("O/C=C\\O", "(Z)-1,2-ethenediol"),
        ("C=C=O", "ethenone"),
        ("C[H](=O)CC", "propanal"),
        ("NN=N", "triazene"),
        ("CCC#C", "but-1-yne"),
    ]
    
    for smi, name in test_cases:
        result, reason = is_gas_molecular_entity(smi)
        print(f"SMILES: {smi} | NAME: {name} -> {result}: {reason}")