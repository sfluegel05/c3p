"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:25438 nucleobase analogue
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for nucleobase cores and analogues
    nucleobase_cores = [
        # Adenine and analogues
        Chem.MolFromSmarts('n1(c)cnc2ncnc12'),          # Adenine core
        Chem.MolFromSmarts('n1cnc2c(n1)[nH]cn2'),       # 2-Aminopurine
        # Guanine and analogues
        Chem.MolFromSmarts('O=C1NC=NC2=NC=NC12'),       # Guanine core
        Chem.MolFromSmarts('O=C1NC(=O)C2=C1N=C(N)N2'),  # Xanthine
        Chem.MolFromSmarts('O=C1NC(=O)C2=C1NC=NC2'),    # Hypoxanthine
        # Cytosine and analogues
        Chem.MolFromSmarts('O=C1NC=CN=C1N'),            # Cytosine core
        Chem.MolFromSmarts('O=C1N=CN=C(C1)N'),          # 5-Substituted cytosine
        # Uracil and analogues
        Chem.MolFromSmarts('O=C1NC(=O)C=CC1'),          # Uracil core
        Chem.MolFromSmarts('O=C1NC(=O)C=C[NH]1'),       # 6-Azauracil
        Chem.MolFromSmarts('O=C1NC(=O)C=CN1C'),         # Thymine core
        # Modified nucleobases
        Chem.MolFromSmarts('n1cnc2c1ncnc2=O'),          # 8-Oxoadenine
        Chem.MolFromSmarts('n1c(=O)[nH]c2c1ncnc2'),     # 8-Hydroxyadenine
        Chem.MolFromSmarts('c1ncnc2[nH]ncnc12'),        # 8-Azaadenine
    ]

    # Check if the molecule matches any of the nucleobase cores
    is_nucleobase = False
    for core in nucleobase_cores:
        if mol.HasSubstructMatch(core):
            is_nucleobase = True
            break

    if not is_nucleobase:
        return False, "Does not contain nucleobase core"

    # Additional checks to filter out large or complex molecules
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_heavy_atoms > 30:
        return False, f"Too many heavy atoms ({num_heavy_atoms}) for a nucleobase analogue"

    num_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 3:
        return False, f"Too many rings ({num_rings}) for a nucleobase analogue"

    mol_weight = Descriptors.ExactMolWt(mol)
    if mol_weight > 400:
        return False, f"Molecular weight ({mol_weight:.2f}) too high for a nucleobase analogue"

    # Check for unusual elements (exclude metals and uncommon non-metals)
    allowed_atomic_nums = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}  # H, C, N, O, F, P, S, Cl, Br, I
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains unusual element ({atom.GetSymbol()}) not typical of nucleobase analogues"

    return True, "Matches nucleobase analogue pattern"