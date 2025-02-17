"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: Gas Molecular Entity
Definition: Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa)
This classifier uses a heuristic based on the allowed elements and a strict limit on the number of heavy atoms.
It also makes two further adjustments:
 - Monoatomic species are accepted only if they are noble gases.
 - If a molecule has no heavy atoms (i.e. only hydrogens) then it is accepted only if it contains exactly two atoms.
 - For molecules with 4 heavy atoms, a combination of S + O together with ≥2 halogen atoms (F, Cl, Br, I)
   is taken as a warning sign that the molecule is likely to be a liquid/solid.
Note: This is just one heuristic approach.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a given molecule (specified as a SMILES string) is a gas molecular entity at STP.
    
    Heuristic criteria used:
      1. The SMILES must represent a valid, single connected molecule.
      2. Every atom must be a main group element (allowed atomic numbers: 1, 2, 5,6,7,8,9,10,11,12,
         13,14,15,16,17,18,19,20,35,36,53,54,86).
      3. If the molecule is monoatomic then it is accepted only if it is a noble gas.
      4. Molecules composed solely of hydrogen isotopes (zero “heavy atoms”) are accepted only if they
         consist of exactly two atoms (e.g. H2, [1H][1H], [3H][3H]).
      5. For a bound molecular entity built from more than one atom, we restrict to very “small” molecules,
         accepting only those with at most 4 heavy atoms.
      6. An additional check for molecules with exactly 4 heavy atoms: if the molecule contains both sulfur (S)
         and oxygen (O) and at least 2 halogen atoms (F, Cl, Br, I) then it is rejected (this catches e.g. thionyl chloride).
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a gas molecular entity, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is a single connected species.
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Molecule is disconnected; should be a single molecular entity"
    
    # Allowed main group atomic numbers: 
    allowed_atomic_numbers = {1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                              14, 15, 16, 17, 18, 19, 20, 35, 36, 53, 54, 86}
    # Noble gases (by symbol)
    noble_gases = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}
    
    # Check that every atom is allowed.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_numbers:
            return False, f"Atom {atom.GetSymbol()} (atomic number {atom.GetAtomicNum()}) is not a main group element"
    
    # Count total atoms and heavy atoms (atomic number > 1)
    total_atoms = mol.GetNumAtoms()
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_heavy = len(heavy_atoms)
    
    # Special case: monoatomic molecules.
    if total_atoms == 1:
        atom = mol.GetAtomWithIdx(0)
        if atom.GetSymbol() in noble_gases:
            mw = Descriptors.ExactMolWt(mol)
            return True, f"Monoatomic noble gas ({atom.GetSymbol()}, MW {mw:.1f} Da) is classified as a gas"
        else:
            return False, f"Monoatomic {atom.GetSymbol()} is not typically a gas at STP"
    
    # Special case: molecules with no heavy atoms (only H isotopes).
    if n_heavy == 0:
        if total_atoms == 2:
            mw = Descriptors.ExactMolWt(mol)
            return True, f"Molecule composed solely of hydrogens (MW {mw:.1f} Da) is classified as a gas"
        else:
            return False, "Molecule composed solely of hydrogen isotopes and with unexpected atom count"
    
    # For multi-atom molecules (n_heavy >= 1), impose a cutoff on heavy atom count.
    if n_heavy > 4:
        return False, f"Too many heavy atoms ({n_heavy}); typical gas molecular entities have ≤4 heavy atoms"
    
    # For cases with exactly 4 heavy atoms, apply an additional check:
    if n_heavy == 4:
        # Count occurrences of S, O, and halogen atoms (F, Cl, Br, I)
        symbols = [atom.GetSymbol() for atom in heavy_atoms]
        count_S = symbols.count("S")
        count_O = symbols.count("O")
        halogens = {"F", "Cl", "Br", "I"}
        count_halogen = sum(1 for sym in symbols if sym in halogens)
        # If we have both S and O and at least 2 halogens, reject.
        if count_S >= 1 and count_O >= 1 and count_halogen >= 2:
            return False, "Molecule contains S, O, and multiple halogens; likely a condensed (non‐gas) entity"
    
    # Otherwise, compute some basic properties for the explanation.
    mw = Descriptors.ExactMolWt(mol)
    # For reporting, we calculate a “halogen fraction” among heavy atoms.
    halogen_atomic_numbers = {9, 17, 35, 53}  # F, Cl, Br, I
    n_halogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() in halogen_atomic_numbers)
    halogen_fraction = n_halogen / n_heavy if n_heavy else 0
    
    return True, (f"Molecule (MW {mw:.1f} Da, {n_heavy} heavy atoms, halogen fraction {halogen_fraction:.2f}) "
                  "is classified as a gas molecular entity")

# Example usage (uncomment the lines below to test):
# test_smiles = [
#     "[H]\\C(C)=C(\\[H])C",  # trans-but-2-ene
#     "FC=C",                # fluoroethene
#     "[H]N([H])[H]",        # ammonia
#     "Cl[H]",               # hydrogen chloride
#     "[H]C([H])([H])[H]",    # methane
#     "[3He]",               # helium-3 atom
#     "C=C",                 # ethene
#     "[1H][1H]",            # diprotium
#     "[220Rn]",             # radon-220 atom
#     "[Rn]",                # radon(0)
#     "[H][H]",              # dihydrogen
#     "[He]",                # helium(0)
#     "FF",                  # difluorine
#     "[6He]",               # helium-6 atom
#     "CC",                  # ethane
#     "ClC=C",               # chloroethene
#     "[219Rn]",             # radon-219 atom
#     "ClCl",                # dichlorine
#     "[Xe]",                # xenon atom
#     "[H]\\C(C)=C(/[H])C",   # cis-but-2-ene
#     "I[H]",                # hydrogen iodide
#     "CCC=C",               # but-1-ene
#     "CC#C",                # propyne
#     "CC(C)C",              # isobutane
#     "[H]C(C)=C([H])C",      # but-2-ene
#     "[O][O]",              # triplet dioxygen
#     "[C-]#[O+]",           # carbon monoxide
#     "[4He]",               # helium-4 atom
#     "[Kr]",                # krypton atom
#     "C1CO1",               # oxirane
#     "[Ar]",                # argon atom
#     "O=C=O",               # carbon dioxide
#     "[222Rn]",             # radon-222 atom
#     "O=[13C]=O",           # ((13)C)carbon dioxide
#     "[3H][3H]",            # ditritium
#     "[C]",                 # carbon atom
# ]
# for smi in test_smiles:
#     result, reason = is_gas_molecular_entity(smi)
#     print(f"SMILES: {smi} -> {result}: {reason}")