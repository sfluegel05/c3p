"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: Gas Molecular Entity
Definition: Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0°C and 100 kPa)
This improved classifier applies tighter criteria based on molecular weight and the number of heavy atoms.
It also adjusts for highly halogenated compounds.
Note: This classifier is heuristic and may not be 100% accurate.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a given molecule (specified as a SMILES string) is a gas molecular entity at STP.
    
    Heuristic criteria used:
      1. The SMILES must represent a valid, connected molecule.
      2. Every atom must be a main group element (allowed: H, He, B, C, N, O, F, Ne, Na, Mg,
         Al, Si, P, S, Cl, Ar, K, Ca, Br, Kr, I, Xe, Rn).
      3. The exact molecular weight must be below 300 Da.
      4. The molecule must be “small” – we use a tight cutoff on the number of heavy atoms 
         (atoms with atomic number > 1). For most molecules we allow at most 7 heavy atoms.
      5. However, if the molecule is heavily halogenated (F, Cl, Br, I), we allow a higher heavy–atom limit.
         Specifically, if more than 50% of the heavy atoms are halogens, we allow up to 15 heavy atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a gas molecular entity, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Allowed main group atomic numbers:
    allowed_atomic_numbers = {1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                              14, 15, 16, 17, 18, 19, 20, 35, 36, 53, 54, 86}
    # Check each atom
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_numbers:
            return False, f"Atom {atom.GetSymbol()} (atomic number {atom.GetAtomicNum()}) is not a main group element"
    
    # Ensure the molecule is a single connected species
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Molecule is disconnected; should be a single molecular entity"

    # Calculate exact molecular weight
    mw = Descriptors.ExactMolWt(mol)
    if mw >= 300:
        return False, f"Molecular weight {mw:.1f} Da is too high for a typical gas molecular entity"
    
    # Count heavy atoms (non-hydrogen atoms) and halogen atoms
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_heavy = len(heavy_atoms)
    
    # Define halogen atomic numbers for F, Cl, Br, I
    halogen_atomic_numbers = {9, 17, 35, 53}
    n_halogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() in halogen_atomic_numbers)
    
    # If there are no heavy atoms then something is off
    if n_heavy == 0:
        return False, "No heavy atoms found; not a typical molecular entity"
    
    # Determine maximum allowed heavy atoms:
    # If more than 50% of heavy atoms are halogens, allow up to 15 heavy atoms (e.g., perfluorocarbons).
    # Otherwise, require at most 7 heavy atoms.
    halogen_fraction = n_halogen / n_heavy
    if halogen_fraction >= 0.5:
        max_heavy_allowed = 15
    else:
        max_heavy_allowed = 7
    
    if n_heavy > max_heavy_allowed:
        return False, f"Too many heavy atoms ({n_heavy}); maximum allowed is {max_heavy_allowed} for a gaseous molecule"
    
    return True, f"Molecule (MW {mw:.1f} Da, {n_heavy} heavy atoms, halogen fraction {halogen_fraction:.2f}) is classified as a gas molecular entity"

# Example usage (uncomment to test):
# examples = [
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
# ]
# for smi in examples:
#     result, reason = is_gas_molecular_entity(smi)
#     print(f"SMILES: {smi} -> {result}: {reason}")