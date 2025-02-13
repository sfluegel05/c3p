"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    An 11,12-saturated fatty acyl-CoA(4-) is a fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts('C(C)(C)(CO)C(=O)NCCC(=O)NCCS.CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA substructure found"

    # Find fatty acyl chain
    acyl_chain = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIsAromatic() and atom.GetDegree() == 3:
            acyl_chain = atom
            break

    if acyl_chain is None:
        return False, "No fatty acyl chain found"

    # Extract fatty acyl chain as a substructure
    acyl_mol = Chem.PathToSubmol(mol, Chem.FindAtomEnvironmentOfRadiusN(mol, acyl_chain.GetIdx(), 10), atomMap=False)

    # Check for 11-12 saturation
    bond_orders = [bond.GetBondType() for bond in acyl_mol.GetBonds()]
    if Chem.BondType.DOUBLE in bond_orders[10:13]:
        return False, "The 11-12 bond of the fatty acyl group is not saturated"

    # Check for fatty acyl chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(acyl_mol)
    if n_rotatable < 6:
        return False, "Fatty acyl chain is too short"

    # Check for presence of double bonds in the fatty acyl chain
    if Chem.BondType.DOUBLE not in bond_orders:
        return False, "Fatty acyl chain must have at least one double bond"

    # Check for stereochemistry at the 3-position
    acyl_atoms = list(acyl_mol.GetAtoms())
    if acyl_atoms[2].GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "Incorrect stereochemistry at the 3-position of the fatty acyl chain"

    return True, "The molecule is an 11,12-saturated fatty acyl-CoA(4-)"