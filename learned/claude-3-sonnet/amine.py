"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule contains an amine group based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two 
    or three hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains amine group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check if molecule contains nitrogen
    if not any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "No nitrogen atoms found"

    # SMARTS patterns for different types of non-amine nitrogen compounds
    exclusion_patterns = {
        "amide": "[NX3]([CX3](=[OX1]))",  # Amide
        "nitro": "[NX3](=[OX1])(=[OX1])",  # Nitro group
        "nitroso": "[NX2]=[OX1]",  # Nitroso
        "azide": "[NX2]=[NX2]=[NX1-,NX2+]",  # Azide
        "diazo": "[NX2]=[NX2]",  # Diazo
        "nitrile": "[NX1]#[CX2]",  # Nitrile
        "isocyanate": "[NX2]=[CX2]=[OX1]",  # Isocyanate
        "isothiocyanate": "[NX2]=[CX2]=[SX1]",  # Isothiocyanate
    }

    # Check for amine patterns
    amine_patterns = {
        "primary": "[NX3H2][CX4]",  # Primary amine (RNH2)
        "secondary": "[NX3H1]([CX4])[CX4]",  # Secondary amine (R2NH)
        "tertiary": "[NX3]([CX4])([CX4])[CX4]",  # Tertiary amine (R3N)
        "aromatic_amine": "[NX3H2]c", # Primary aromatic amine
        "aromatic_secondary": "[NX3H1](c)[CX4]", # Secondary with one aromatic
        "aromatic_tertiary": "[NX3](c)([CX4])[CX4]" # Tertiary with aromatic
    }

    # Check if any nitrogen is part of excluded groups
    for name, pattern in exclusion_patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(pat):
            # Don't immediately return False - the molecule might have other valid amine groups
            continue

    # Check for amine patterns
    found_amine = False
    amine_types = []
    
    for name, pattern in amine_patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(pat):
            found_amine = True
            amine_types.append(name)

    if not found_amine:
        return False, "No amine groups found"

    # Additional validation: Check that the nitrogen is not part of a ring system
    # that would make it not a true amine (like pyridine)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            if atom.GetIsAromatic():  # Skip aromatic nitrogens in rings
                continue
            # Check if it's bonded to carbon with single bonds
            carbon_single_bonds = sum(1 for bond in atom.GetBonds() 
                                   if bond.GetBondType() == Chem.BondType.SINGLE 
                                   and bond.GetOtherAtom(atom).GetAtomicNum() == 6)
            if carbon_single_bonds > 0:
                found_amine = True

    if not found_amine:
        return False, "Contains nitrogen but not as amine"

    amine_description = ", ".join(amine_types)
    return True, f"Contains amine groups of types: {amine_description}"