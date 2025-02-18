"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Minerals including ionic compounds, covalent structures, and elemental forms formed through geological processes.
"""
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Valid minerals include:
    - Ionic compounds with metal cations and common inorganic anions
    - Covalent metal oxides/sulfides/halides
    - Elemental metals/metalloids or sulfur/carbon (graphite/diamond)
    - Hydrates of the above

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if mineral, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define metals/metalloids and allowed nonmetals for covalent structures
    metals = {'Li','Be','Na','Mg','Al','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',
              'Ni','Cu','Zn','Ga','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
              'Ag','Cd','In','Sn','Sb','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
              'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir',
              'Pt','Au','Hg','Tl','Pb','Bi','Po','Fr','Ra','Ac','Th','Pa','U','Si','B','As','Se','Te','Ge'}
    
    # Common inorganic anion patterns (carbonate, sulfate, phosphate, etc.)
    inorganic_anion_patterns = [
        Chem.MolFromSmarts('[O-]S(=O)(=O)[O-]'),  # Sulfate
        Chem.MolFromSmarts('[O-]P(=O)([O-])[O-]'), # Phosphate
        Chem.MolFromSmarts('[O-]C(=O)[O-]'),       # Carbonate
        Chem.MolFromSmarts('[O-][N+](=O)[O-]'),    # Nitrate
        Chem.MolFromSmarts('[O-]Cl'),              # Hypochlorite
        Chem.MolFromSmarts('[S-]'),                # Sulfide
        Chem.MolFromSmarts('[F-]'),                # Fluoride
        Chem.MolFromSmarts('[Cl-]'),               # Chloride
        Chem.MolFromSmarts('[OH-]'),               # Hydroxide
        Chem.MolFromSmarts('[O-][Si]'),            # Silicate groups
        Chem.MolFromSmarts('[O-][B]')              # Borate groups
    ]

    # Split into fragments (ions/metal complexes)
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    
    # Remove water fragments (hydrates allowed)
    non_water_frags = []
    for frag in fragments:
        formula = rdMolDescriptors.CalcMolFormula(frag)
        if formula == 'H2O':
            continue
        non_water_frags.append(frag)
    
    # Check 1: Ionic compounds (multiple fragments with metal cations and inorganic anions)
    if len(non_water_frags) >= 2:
        has_metal_cation = False
        has_inorganic_anion = False
        
        for frag in non_water_frags:
            frag_charge = sum(atom.GetFormalCharge() for atom in frag.GetAtoms())
            contains_metal = any(atom.GetSymbol() in metals for atom in frag.GetAtoms())
            contains_carbon = any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms())
            
            # Check for metal cation
            if frag_charge > 0 and contains_metal:
                has_metal_cation = True
            
            # Check for inorganic anion (allow carbon only in specific groups)
            if frag_charge < 0:
                if any(frag.HasSubstructMatch(patt) for patt in inorganic_anion_patterns):
                    has_inorganic_anion = True
                elif not contains_carbon:  # Allow non-carbon anions
                    has_inorganic_anion = True
        
        if has_metal_cation and has_inorganic_anion:
            return True, "Ionic compound with metal cation and inorganic anion"

    # Check 2: Covalent structures (single fragment with metal + nonmetal bonds)
    if len(non_water_frags) == 1:
        has_metal = False
        has_mineral_nonmetal = False
        for atom in non_water_frags[0].GetAtoms():
            sym = atom.GetSymbol()
            if sym in metals:
                has_metal = True
            elif sym in {'O','S','F','Cl','Br','I'}:
                has_mineral_nonmetal = True
        
        # Require direct metal-nonmetal bonds
        if has_metal and has_mineral_nonmetal:
            # Check for at least one metal-nonmetal bond
            for bond in non_water_frags[0].GetBonds():
                a1 = bond.GetBeginAtom().GetSymbol()
                a2 = bond.GetEndAtom().GetSymbol()
                if (a1 in metals and a2 in {'O','S','F','Cl','Br','I'}) or \
                   (a2 in metals and a1 in {'O','S','F','Cl','Br','I'}):
                    return True, "Covalent metal-nonmetal compound"

    # Check 3: Elemental forms (metals/metalloids or sulfur/carbon)
    elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if len(elements) == 1:
        elem = elements.pop()
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        if total_charge != 0:
            return False, "Charged elemental form"
        if elem in metals or elem in {'S','C'}:
            return True, "Elemental form of a mineral"

    return False, "Does not meet mineral criteria"