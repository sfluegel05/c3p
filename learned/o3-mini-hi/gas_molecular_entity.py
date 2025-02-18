"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: Gas Molecular Entity
Definition: Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa).

This heuristic classifier works as follows:
  1. The SMILES is parsed and explicit hydrogens are added so that we count atoms properly.
  2. For monoatomic molecules (after adding hydrogens) only noble gases (and their isotopes) are accepted.
  3. For polyatomic molecules:
       – Every heavy atom (atomic number > 1) must be in the safe set: {C, N, O, F, Cl, Br, I}.
       – The molecule must be neutral (all atoms with formal charge 0).
       – The number of heavy atoms must be small (≤4).
       – Extra heuristic rejections:
             * If the number of heavy atoms is 0 (only hydrogens) then exactly 2 atoms are allowed.
             * If (after adding Hs) the heavy–atom count is 2 and the set equals {"C","N"} we reject (e.g. methanimine).
             * If the heavy–atom count is 4 and the set is exactly {"C", "N"} (as in trimethylamine) we reject.
             * If the heavy–atom count is ≥3 and the fraction of heavy atoms that are halogens (F,Cl,Br,I) is 0.67 or greater then reject
               (except if it is a diatomic halogen where both atoms are the same, e.g. FF or ClCl).
             * If any triple bond exists between a carbon and a nitrogen (a nitrile motif) then reject.
             * If the molecule has the isocyanate motif “[N]=C=O” then reject.
  4. Otherwise the molecule is accepted as a gas molecular entity.
    
Note: This heuristic uses structure‐patterns as surrogates for likely boiling points.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a given molecule (specified as a SMILES string) is a gas molecular entity at STP
    based on a series of structural heuristics.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a gas molecular entity, False otherwise.
        str: Explanation/reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that implicit H atoms are considered in counting.
    mol = Chem.AddHs(mol)
    
    # Ensure the molecule is a single connected species.
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Molecule is disconnected; should be a single connected entity"

    total_atoms = mol.GetNumAtoms()  # including explicit hydrogens
    
    # Separate heavy atoms (atomic number > 1) from hydrogens.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_heavy = len(heavy_atoms)
    
    # Special case: monoatomic species.
    # Note: For monatomic species, implicit Hs are not added.
    if total_atoms == 1:
        atom = mol.GetAtomWithIdx(0)
        symbol = atom.GetSymbol()
        noble_gases = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}
        if symbol in noble_gases:
            mw = Descriptors.ExactMolWt(mol)
            return True, f"Monoatomic noble gas ({symbol}, MW {mw:.1f} Da) is classified as a gas"
        else:
            return False, f"Monoatomic {symbol} is not classified as a gas (only noble gases accepted as atoms)"
    
    # Special case: molecules composed solely of hydrogens
    if n_heavy == 0:
        if total_atoms == 2:
            mw = Descriptors.ExactMolWt(mol)
            return True, f"Molecule composed solely of hydrogens (MW {mw:.1f} Da) is classified as a gas"
        else:
            return False, "Molecule composed solely of hydrogen isotopes but with unexpected atom count"
    
    # For polyatomic species, check that every heavy atom is from the safe set.
    safe_heavy = {"C", "N", "O", "F", "Cl", "Br", "I"}
    heavy_symbols = []
    for atom in heavy_atoms:
        if atom.GetFormalCharge() != 0:
            return False, f"Atom {atom.GetSymbol()} carries nonzero formal charge; molecule rejected"
        sym = atom.GetSymbol()
        heavy_symbols.append(sym)
        if sym not in safe_heavy:
            return False, f"Atom {sym} is not allowed for gas molecular entities"
    
    # Enforce maximum number of heavy atoms.
    if n_heavy > 4:
        return False, f"Too many heavy atoms ({n_heavy}); typical gas molecular entities have ≤4 heavy atoms"
    
    hs_set = set(heavy_symbols)
    
    # Ad hoc rejection rules based on heavy–atom composition:
    if n_heavy == 2 and hs_set == {"C", "N"}:
        return False, f"Molecule has 2 heavy atoms with composition {hs_set} (e.g. methanimine); likely not a gas"
    if n_heavy == 4 and hs_set == {"C", "N"}:
        return False, "Molecule has 4 heavy atoms composed solely of C and N (e.g. trimethylamine); likely not a gas"
    
    # Compute crude halogen fraction (F, Cl, Br, I) among heavy atoms.
    halogen_atomic_nums = {9, 17, 35, 53}
    n_halogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() in halogen_atomic_nums)
    halogen_fraction = n_halogen / n_heavy if n_heavy else 0

    # If there are 3 or more heavy atoms and a very high halogen fraction (>=0.67), reject—
    # except if the molecule is a diatomic halogen of one element (which are accepted).
    if n_heavy >= 3 and halogen_fraction >= 0.67:
        return False, f"High halogen fraction ({halogen_fraction:.2f}) with {n_heavy} heavy atoms; likely not a gas"
    
    # For diatomic molecules composed solely of halogens, ensure they are homodiatomic.
    if n_heavy == 2 and all(atom.GetAtomicNum() in halogen_atomic_nums for atom in heavy_atoms):
        if heavy_symbols[0] != heavy_symbols[1]:
            return False, f"Heteronuclear diatomic halogen ({heavy_symbols[0]},{heavy_symbols[1]}) not accepted as a gas"
    
    # Check for a triple bond between carbon and nitrogen (nitrile motif).
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            syms = {a1.GetSymbol(), a2.GetSymbol()}
            if "C" in syms and "N" in syms:
                return False, "Molecule has a C≡N triple bond (nitrile motif); likely not a gas"
    
    # Check for isocyanate functional group [N]=C=O.
    pat_iso = Chem.MolFromSmarts("[N]=C=O")
    if pat_iso and mol.HasSubstructMatch(pat_iso):
        return False, "Molecule contains the isocyanate motif [N]=C=O; likely not a gas"
    
    # If we passed all tests, accept the molecule.
    mw = Descriptors.ExactMolWt(mol)
    return True, (f"Molecule (MW {mw:.1f} Da, {n_heavy} heavy atoms, halogen fraction {halogen_fraction:.2f}) "
                  "is classified as a gas molecular entity")

# Example usage (you can remove or modify this section when integrating into a larger project):
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
        # False negatives (must be accepted but were missed before)
        ("[H]N([H])[H]", "ammonia"),
        ("Cl[H]", "hydrogen chloride"),
        ("[H]C([H])([H])[H]", "methane"),
        ("I[H]", "hydrogen iodide"),
        # Some borderline cases also
        ("[C-]#[O+]", "carbon monoxide"),
        ("[C]", "carbon atom"),
        ("[O-][O+]=O", "ozone"),
        ("C1OC=C1", "oxetene"),
        ("FCCl", "Chlorofluoromethane"),
    ]
    
    for smi, name in test_cases:
        result, reason = is_gas_molecular_entity(smi)
        print(f"SMILES: {smi} NAME: {name} -> {result}: {reason}")