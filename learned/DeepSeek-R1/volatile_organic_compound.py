"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: Volatile Organic Compound (VOC) based on boiling point estimation
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on SMILES.
    VOC definition: Organic compounds with boiling point <=250°C at 101.3 kPa.
    Uses molecular weight and functional groups as boiling point proxies.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Must contain carbon to be organic
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms"

    # Calculate molecular properties
    mol_wt = Descriptors.ExactMolWt(mol)
    h_bond_donors = Chem.rdMolDescriptors.CalcNumHBD(mol)
    rotatable_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Main chain length estimation (crude approximation)
    longest_carbon_chain = max(
        (len(match) for match in mol.GetSubstructMatches(Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]"))),
        default=0
    )

    # Empirical rules for boiling point estimation
    # Rule 1: Short molecules are always VOCs
    if mol_wt <= 150:
        return True, f"Molecular weight {mol_wt:.1f} <=150"

    # Rule 2: Long carbon chains (>12 carbons) typically exceed 250°C
    if longest_carbon_chain > 12:
        return False, f"Long chain ({longest_carbon_chain} carbons)"

    # Rule 3: Heavy molecules (>300 g/mol) are excluded
    if mol_wt > 300:
        return False, f"Molecular weight {mol_wt:.1f} >300"

    # Rule 4: Molecules with H-bond donors get stricter checks
    if h_bond_donors > 0:
        # Alcohols with >8 carbons typically exceed 250°C
        if longest_carbon_chain > 8:
            return False, f"Long chain alcohol ({longest_carbon_chain} carbons)"
        # Heavy polar molecules (>200 g/mol) excluded
        if mol_wt > 200:
            return False, f"Polar molecule weight {mol_wt:.1f} >200"

    # Rule 5: Highly branched molecules get more leniency
    if rotatable_bonds < 3 and mol_wt <= 250:
        return True, "Branched structure with moderate weight"

    # Fallback: Accept if passes basic checks
    return True, "Passed empirical volatility checks"