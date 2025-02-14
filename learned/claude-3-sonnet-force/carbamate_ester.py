"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:35669 carbamate ester
'Any ester of carbamic acid or its N-substituted derivatives.'
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbamate functional group (-O-C(=O)-N-)
    carbamate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[NX3]")
    carbamate_matches = mol.GetSubstructMatches(carbamate_pattern)

    if not carbamate_matches:
        return False, "No carbamate functional group found"

    # Check for ester linkage on nitrogen
    ester_pattern = Chem.MolFromSmarts("[NX3][OX2][CX3]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "Carbamate nitrogen not esterified"

    # Check for allowed substituents on nitrogen
    allowed_substituents = ["C", "c", "O", "N", "S", "P", "F", "Cl", "Br", "I"]
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    for n_atom in n_atoms:
        substituents = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in [mol.GetBondWithAtomAsIndex(n_atom.GetIdx(), i).GetOtherAtomIdx(n_atom.GetIdx()) for i in range(n_atom.GetNumExplicitNeighbors())]]
        if any(sub not in allowed_substituents for sub in substituents):
            return False, "Disallowed substituent on carbamate nitrogen"

    # Exclude molecules too small to be carbamate esters
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for carbamate ester"

    return True, "Contains carbamate functional group with ester linkage and allowed substituents"