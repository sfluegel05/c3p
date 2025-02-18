"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Mineral
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.

    Minerals are naturally occurring substances formed through geological processes.
    This function uses heuristics to classify minerals based on the absence of organic
    functional groups, the presence of elements commonly found in minerals, and the
    presence of inorganic anions and cations.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is likely a mineral, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # List of elements commonly found in minerals (atomic numbers)
    mineral_elements = set([
        1,  # Hydrogen
        3, 4, 11, 12, 13, 14,  # Li, Be, Na, Mg, Al, Si
        15, 16, 17, 19, 20, 25, 26,  # P, S, Cl, K, Ca, Mn, Fe
        27, 28, 29, 30,  # Co, Ni, Cu, Zn
        33, 34, 35,  # As, Se, Br
        37, 38, 39,  # Rb, Sr, Y
        47, 48,  # Ag, Cd
        53, 54, 55, 56,  # I, Xe, Cs, Ba
        79, 80, 81, 82, 83,  # Au, Hg, Tl, Pb, Bi
        92,  # Uranium
        # Add more elements as needed
    ])

    # Check for presence of elements commonly found in minerals
    elements_in_mol = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    has_mineral_element = any(elem in mineral_elements for elem in elements_in_mol)
    if not has_mineral_element:
        return False, "No elements commonly found in minerals detected"

    # Check for presence of organic functional groups
    organic_patterns = [
        '[CX4]',  # sp3 Carbon
        '[c]',    # Aromatic carbon
        '[#6]=[#6]',  # Carbon double bond
        '[#6]#[#6]',  # Carbon triple bond
        '[#6]-[#1]',  # Carbon-Hydrogen bond
        '[#7]-[#1]',  # Nitrogen-Hydrogen bond
        '[#8]-[#1]',  # Oxygen-Hydrogen bond
    ]
    for pattern in organic_patterns:
        smarts = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(smarts):
            return False, "Contains organic functional groups"

    # Exclude molecules with long carbon chains (indicative of organic compounds)
    num_chains = Chem.rdmolops.GetSSSR(mol)
    if num_chains > 0:
        longest_chain = max(len(chain) for chain in Chem.rdmolops.GetRingInfo(mol).AtomRings())
        if longest_chain > 6:
            return False, "Contains long carbon chains"

    # Include common inorganic anions and cations in minerals
    inorganic_ions = [
        '[O-]C([O-])=O',             # Carbonate
        '[O-]S(=O)(=O)[O-]',         # Sulfate
        '[S-2]',                     # Sulfide
        '[N-3]',                     # Nitride
        '[O-]', '[OH-]',             # Oxide, Hydroxide
        '[Cl-]', '[Br-]', '[I-]',    # Halides
        'P(=O)([O-])([O-])[O-]',     # Phosphate
        '[F-]',                      # Fluoride
        '[Si]', '[SiH4]',            # Silicates (simplified)
        '[N+](=O)[O-]',              # Nitrate
        '[C-]#[N+]',                 # Cyanide
        '[As]', '[AsH3]',            # Arsenic compounds (simplified)
    ]
    has_inorganic_ion = False
    for ion in inorganic_ions:
        ion_pattern = Chem.MolFromSmarts(ion)
        if mol.HasSubstructMatch(ion_pattern):
            has_inorganic_ion = True
            break
    if not has_inorganic_ion:
        return False, "No common inorganic ions detected"

    # Check for absence of complex organic structures
    if rdMolDescriptors.CalcNumAromaticRings(mol) > 0:
        return False, "Contains aromatic rings (indicative of organic compounds)"

    # Check for presence of metal-nonmetal bonds
    metals = set([
        3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27,
        28, 29, 30, 31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46,
        47, 48, 49, 55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79,
        80, 81, 82, 83, 84, 87, 88, 89, 104, 105, 106, 107, 108,
        109, 110, 111, 112
    ])
    nonmetals = set([
        1, 6, 7, 8, 9, 15, 16, 17, 34, 35, 53,  # H, C, N, O, F, P, S, Cl, Se, Br, I
    ])
    has_metal_nonmetal_bond = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if ((a1.GetAtomicNum() in metals and a2.GetAtomicNum() in nonmetals) or
            (a2.GetAtomicNum() in metals and a1.GetAtomicNum() in nonmetals)):
            has_metal_nonmetal_bond = True
            break
    if not has_metal_nonmetal_bond:
        return False, "No metal-nonmetal bonds detected"

    return True, "Likely a mineral (inorganic compound with elements commonly found in minerals)"