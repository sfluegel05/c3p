"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base
    in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to ensure aromaticity is perceived correctly
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass

    # Identify nucleobase (purine or pyrimidine ring system)
    # More general patterns that allow for substitutions
    purine_base = Chem.MolFromSmarts('c1ncnc2n[cH0,n]c[nH0,n]c12')  # purine core
    pyrimidine_base = Chem.MolFromSmarts('c1[nH0,n]c[nH0,n]c[nH0,n]c1')  # pyrimidine core

    has_purine = mol.HasSubstructMatch(purine_base)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_base)

    if not (has_purine or has_pyrimidine):
        return False, "Nucleobase not found (purine or pyrimidine base)"

    # Identify sugar (ribose or deoxyribose with possible modifications)
    # General pattern for a furanose ring (five-membered ring with oxygen)
    sugar = Chem.MolFromSmarts('C1OC([#6,H])C([#6,H])C1')  # furanose ring

    has_sugar = mol.HasSubstructMatch(sugar)
    if not has_sugar:
        return False, "Sugar (ribose or deoxyribose) not found"

    # Check for N-glycosidic bond between nucleobase and sugar
    # Find an N-C bond between the base and the sugar
    nucleobase_nitrogens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.IsAromatic()]
    sugar_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.IsInRingSize(5)]

    glycosidic_bond_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if ((a1 in nucleobase_nitrogens and a2 in sugar_carbons) or
            (a2 in nucleobase_nitrogens and a1 in sugar_carbons)):
            if bond.GetBondType() == rdchem.BondType.SINGLE:
                glycosidic_bond_found = True
                break

    if not glycosidic_bond_found:
        return False, "N-glycosidic bond between nucleobase and sugar not found"

    # Check for phosphate group(s) attached to the 5' carbon of the sugar
    # Phosphate group SMARTS pattern
    phosphate = Chem.MolFromSmarts('OP(=O)(O)O')  # phosphate group

    # 5' carbon is bonded to oxygen which is bonded to phosphate
    phosphate_attachment = Chem.MolFromSmarts('C(O[P](=O)(O)O)')  # primary alcohol with phosphate

    # Check for attachment to the sugar ring
    phosphate_bonded = False
    for match in mol.GetSubstructMatches(phosphate_attachment):
        # Check if the carbon is connected to the sugar ring
        phosphate_carbon = mol.GetAtomWithIdx(match[0])
        if any([bond.GetBondType() == rdchem.BondType.SINGLE and
                bond.GetOtherAtomIdx(phosphate_carbon.GetIdx()) in sugar_carbons
                for bond in phosphate_carbon.GetBonds()]):
            phosphate_bonded = True
            break

    if not phosphate_bonded:
        return False, "Phosphate group at 5' position not found"

    # Count the number of phosphate groups attached to the molecule
    phosphate_groups = mol.GetSubstructMatches(phosphate)
    num_phosphates = len(phosphate_groups)
    if num_phosphates < 1 or num_phosphates > 4:
        return False, f"Found {num_phosphates} phosphate groups, but need between 1 and 4"

    return True, "Molecule is a nucleoside 5'-phosphate"

__metadata__ = {
    'chemical_class': {
        'name': "nucleoside 5'-phosphate",
        'definition': "A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated."
    }
}