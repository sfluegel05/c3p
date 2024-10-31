from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkatriene(smiles: str):
    """
    Determines if a molecule is an alkatriene (acyclic hydrocarbon with 3 C=C double bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkatriene, False otherwise
        str: Reason for classification
    """
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check if molecule has rings
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Molecule contains rings - alkatrienes must be acyclic"

    # Count number of double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1
            # Check if double bond is between carbons
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() != 'C' or end_atom.GetSymbol() != 'C':
                return False, "Contains non-carbon-carbon double bonds"

    if double_bond_count != 3:
        return False, f"Contains {double_bond_count} double bonds - alkatrienes must have exactly 3"

    # Check if molecule contains only carbon and hydrogen
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'H']:
            return False, "Contains non-hydrocarbon atoms"

    return True, "Acyclic hydrocarbon with exactly 3 carbon-carbon double bonds"
# Pr=1.0
# Recall=0.96