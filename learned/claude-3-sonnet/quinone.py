"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:38835 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is a compound with a fully conjugated cyclic dione structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enumerate tautomers
    try:
        tautomers = list(enumerate_tautomers(mol))
    except Exception as e:
        return False, f"Failed to enumerate tautomers: {e}"

    # Check for conjugated cyclic dione pattern
    quinone_pattern = Chem.MolFromSmarts("[$(O=C1C=CC=CC=C1=O)]")
    matches = []
    for tautomer in tautomers:
        matches.extend(tautomer.GetSubstructMatches(quinone_pattern))

    if not matches:
        return False, "No conjugated cyclic dione structure found"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for quinones"

    # Check atom counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 6 or o_count < 2:
        return False, "Insufficient carbon or oxygen atoms for quinone"

    return True, "Contains a conjugated cyclic dione structure"

def enumerate_tautomers(mol):
    """
    Generate tautomers for a given molecule.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule

    Yields:
        rdkit.Chem.rdchem.Mol: Tautomer molecules
    """
    tautomer_smiles = [Chem.MolToSmiles(mol)]
    tautomer_mols = [mol]

    while True:
        tautomers = tautomer_mols
        tautomer_mols = []
        for tautomer in tautomers:
            tautomer_mols.extend(AllChem.EnumerateIsomers(tautomer))

        new_tautomers = [Chem.MolToSmiles(mol) for mol in tautomer_mols]
        if not any(smiles not in tautomer_smiles for smiles in new_tautomers):
            break

        tautomer_smiles.extend(new_tautomers)
        tautomer_mols = list(set(tautomer_mols))

    yield from tautomer_mols