"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:35723 carbamate ester
A carbamate ester is an ester derived from carbamic acid (NH2COOH) by
substitution of one or both of the hydrogen atoms of the amino group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Check for carbamate ester pattern: -O-C(=O)-N-
    carbamate_pattern = Chem.MolFromSmarts("[OX2]-C(=O)-[NX3]")
    if not mol.HasSubstructMatch(carbamate_pattern):
        return False, "No carbamate ester substructure found"

    # Check for carbamic acid derivatives (exclude carbazic acid and dimethylcarbamic acid)
    excluded_patterns = [Chem.MolFromSmarts("O=C(N)N"), Chem.MolFromSmarts("CN(C)C(=O)O")]
    for pattern in excluded_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Not a carbamate ester, excluded pattern found"

    # Check for common substituents on nitrogen atom
    allowed_substituents = [Chem.MolFromSmarts("[NH2]"), Chem.MolFromSmarts("[NH1]"), Chem.MolFromSmarts("[NX3H0]")]
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    for n_atom in n_atoms:
        is_allowed = False
        for pattern in allowed_substituents:
            if mol.HasSubstructMatch(pattern, atomIdxs=[n_atom.GetIdx()]):
                is_allowed = True
                break
        if not is_allowed:
            return False, "Nitrogen substituent not allowed for carbamate esters"

    # Check for common molecular scaffolds
    scaffold_patterns = [Chem.MolFromSmarts("c1ccccc1"), Chem.MolFromSmarts("C1CCNCC1"), Chem.MolFromSmarts("C1CCNC1")]
    has_scaffold = any(mol.HasSubstructMatch(pattern) for pattern in scaffold_patterns)

    # Check molecular weight and other descriptors
    mol_wt = AllChem.CalcExactMolWt(mol)
    n_rotatable = AllChem.CalcNumRotatableBonds(mol)
    if mol_wt < 100 or n_rotatable < 2:
        return False, "Molecular weight or rotatable bond count too low for carbamate ester"

    if has_scaffold:
        return True, "Contains carbamate ester substructure and common molecular scaffold"
    else:
        return True, "Contains carbamate ester substructure"