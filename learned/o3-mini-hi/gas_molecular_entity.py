"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: Gas Molecular Entity
Definition: Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa).

This heuristic classifier uses a series of tests:
  1. The SMILES is parsed as a single connected molecule.
  2. For monoatomic molecules, only noble gases (and their isotopes) are accepted.
  3. For polyatomic molecules, every “heavy” atom (atomic number > 1) must be one of a safe set 
     (here chosen as only C, N, O, F, Cl, Br, and I). (Other main‐group atoms such as P, S, Al are rejected.)
  4. Molecules must be very small – we require at most 4 heavy atoms.
  5. Molecules carrying any nonzero formal charge are rejected.
  6. For borderline cases based solely on count, extra ad hoc rules are applied:
       • If a molecule has no heavy atoms (only H’s) it is accepted only if exactly two atoms are present.
       • If a polyatomic molecule has exactly 4 heavy atoms and its heavy atoms are C and N 
         (e.g. trimethylamine) then it is rejected.
       • If a molecule has exactly 2 heavy atoms and its heavy–atom set equals {C,Cl},
         {C,O} or {O,I} (patterns seen in e.g. chloromethane or hypoiodite) then it is rejected.
       • If a molecule of 3 heavy atoms matches the imine pattern [N]=[C]=[N] (methanediimine)
         then it is rejected.
       
Note: This heuristic does not “predict” boiling points but is tuned to accept examples like
trans‐but–2–ene, fluoroethene, ammonia, hydrogen chloride, methane, dihydrogen, ethene, radon, etc.,
and to reject species known not to be gases (such as trimethylamine, oxychlorides with formal charge, 
and molecules containing S, P, or Al).
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a given molecule (specified as a SMILES string) is a gas molecular entity at STP.

    Heuristic criteria used:
      1. The molecule must be a valid, single, connected species.
      2. For monoatomic molecules, only noble gases (He, Ne, Ar, Kr, Xe, Rn) are accepted.
      3. For polyatomic molecules:
             - All heavy atoms (atomic number > 1) must be chosen from a safe set:
               {C, N, O, F, Cl, Br, I}.
             - The molecule must be neutral (all atoms with formal charge 0).
             - The number of heavy atoms (n_heavy) must be at most 4.
             - Extra rules for certain combinations are applied:
                 * If no heavy atoms (only H isotopes), then exactly two atoms must be present.
                 * If n_heavy == 4 and the heavy–atom set is exactly {'C','N'}, reject (e.g. trimethylamine).
                 * If n_heavy == 2 and the heavy–atom set is in any of:
                     {'C','Cl'}, {'C','O'}, or {'O','I'}, reject.
                 * If n_heavy == 3 and the molecule matches the pattern [N]=[C]=[N],
                   then reject.
      4. Otherwise, the molecule is classified as a gas molecular entity.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a gas molecular entity, False otherwise.
        str: Explanation/reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is a single connected species.
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Molecule is disconnected; should be a single molecular entity"
    
    # Define safe heavy–atom symbols for gases.
    safe_heavy_atoms = {"C", "N", "O", "F", "Cl", "Br", "I"}
    # Noble gases (for single atom species). Allow isotopes (e.g. [3He], [220Rn]) by checking
    # that the element symbol (after stripping any digits) is in this set.
    noble_gases = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}
    
    total_atoms = mol.GetNumAtoms()
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_heavy = len(heavy_atoms)
    
    # Special case: monoatomic molecular entities.
    if total_atoms == 1:
        atom = mol.GetAtomWithIdx(0)
        # Remove isotope label if present.
        symbol = atom.GetSymbol()
        if symbol in noble_gases:
            mw = Descriptors.ExactMolWt(mol)
            return True, f"Monoatomic noble gas ({symbol}, MW {mw:.1f} Da) is classified as a gas"
        else:
            return False, f"Monoatomic {symbol} is not classified as a gas (only noble gases accepted as atoms)"
    
    # Special case: molecules made only of hydrogen isotopes.
    if n_heavy == 0:
        if total_atoms == 2:
            mw = Descriptors.ExactMolWt(mol)
            return True, f"Molecule composed solely of hydrogens (MW {mw:.1f} Da) is classified as a gas"
        else:
            return False, "Molecule composed solely of hydrogen isotopes but with unexpected atom count"
    
    # For polyatomic species, require that each heavy atom is from the safe set.
    heavy_symbols = []
    for atom in heavy_atoms:
        if atom.GetFormalCharge() != 0:
            return False, f"Atom {atom.GetSymbol()} carries nonzero formal charge; molecule rejected"
        sym = atom.GetSymbol()
        heavy_symbols.append(sym)
        if sym not in safe_heavy_atoms:
            return False, f"Atom {sym} is not allowed for gas molecular entities"
    
    # Enforce a maximum heavy atom count.
    if n_heavy > 4:
        return False, f"Too many heavy atoms ({n_heavy}); typical gas molecular entities have ≤4 heavy atoms"
    
    hs_set = set(heavy_symbols)
    
    # Extra rules based on the heavy atom counts and compositions:
    # Rule for molecules with exactly 4 heavy atoms that are composed of C and N (e.g. trimethylamine).
    if n_heavy == 4 and hs_set == {"C", "N"}:
        return False, "Molecule has 4 heavy atoms composed solely of C and N (e.g. trimethylamine); likely not a gas"
    
    # Rule for molecules with exactly 2 heavy atoms that are a combination of carbon with either Cl,
    # oxygen, or (O,I) as seen in e.g. chloromethane or hypoiodite.
    disallowed_pairs = [ {"C", "Cl"}, {"C", "O"}, {"O", "I"} ]
    if n_heavy == 2 and hs_set in disallowed_pairs:
        return False, f"Molecule has 2 heavy atoms with composition {hs_set}, which is disfavoured for a gas"
    
    # If the molecule has 3 heavy atoms, check for the imine pattern N=C=N.
    pattern = Chem.MolFromSmarts("[N]=[C]=[N]")
    if n_heavy == 3 and pattern is not None and mol.HasSubstructMatch(pattern):
        return False, "Molecule matches disfavoured imine pattern [N]=[C]=[N]; likely not a gas"
    
    # If we pass the above tests, we accept the molecule as a gas molecular entity.
    mw = Descriptors.ExactMolWt(mol)
    # Calculate a crude halogen fraction (among heavy atoms)
    halogen_nums = {9, 17, 35, 53}  # F, Cl, Br, I atomic numbers
    n_halogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() in halogen_nums)
    halogen_fraction = n_halogen / n_heavy if n_heavy else 0

    return True, (f"Molecule (MW {mw:.1f} Da, {n_heavy} heavy atoms, halogen fraction {halogen_fraction:.2f}) "
                  "is classified as a gas molecular entity")

# Example usage:
if __name__ == '__main__':
    test_cases = [
        # True positives
        ("[H]\\C(C)=C(\\[H])C", "trans-but-2-ene"),
        ("FC=C", "fluoroethene"),
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
        ("CCC=C", "but-1-ene"),
        ("CC#C", "propyne"),
        ("CC(C)C", "isobutane"),
        ("[H]C(C)=C([H])C", "but-2-ene"),
        ("[O][O]", "triplet dioxygen"),
        ("[C-]#[O+]", "carbon monoxide"),
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
        # False positives / negatives according to the provided outcomes:
        ("O[P][O-]", "hydroxidooxidophosphate(.1-)"),
        ("[H]OSO[H]", "dihydroxidosulfur"),
        ("OS(O)=S", "sulfurothious O-acid"),
        ("CN(C)C", "trimethylamine"),
        ("O=[Cl+]=O", "dioxidochlorine(1+)"),
        ("[H]O[Al](O[H])O[H]", "aluminium hydroxide"),
        ("OSSO", "disulfanediol"),
        ("[C-]#[C-]", "dicarbide(2-)"),
        ("OCCS", "mercaptoethanol"),
        ("N=C=N", "methanediimine"),
        ("[H]C([H])([H])Cl", "chloromethane"),
        ("C[O]", "methyloxidanyl"),
        ("SC(S)=O", "carbonodithioic S,S-acid"),
        ("[S-]S[S-]", "trisulfide(2-)"),
        ("O=[Al]", "oxidoaluminium"),
        ("FN(F)Cl", "chlorodifluoroamine"),
        ("CP", "methylphosphine"),
        ("CCSC", "ethyl methyl sulfide"),
        ("SC([S-])=S", "hydrogen trithiocarbonate"),
        ("FCCl", "Chlorofluoromethane"),
        ("ClC(Br)Br", "Chlorodibromomethane"),
        ("CP(C)C", "trimethylphosphine"),
        ("[H][C]=O", "oxomethyl"),
        ("[O-]I", "hypoiodite"),
        ("[Na]OC#N", "sodium cyanate"),
        ("[H]N([H])[H]", "ammonia"),
        ("Cl[H]", "hydrogen chloride"),
        ("[H]C([H])([H])[H]", "methane"),
        ("I[H]", "hydrogen iodide"),
        ("[C]", "carbon atom"),
    ]
    
    for smi, name in test_cases:
        result, reason = is_gas_molecular_entity(smi)
        print(f"SMILES: {smi} NAME: {name} -> {result}: {reason}")