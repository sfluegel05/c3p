"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are endogenous ligands that activate cannabinoid receptors,
    typically composed of long-chain polyunsaturated fatty acids linked to ethanolamine
    or glycerol via amide or ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify amide or ester bonds
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    has_amide = mol.HasSubstructMatch(amide_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)

    if not (has_amide or has_ester):
        return False, "No amide or ester bond found"

    # Identify head groups
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    has_ethanolamine = mol.HasSubstructMatch(ethanolamine_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)

    if not (has_ethanolamine or has_glycerol):
        return False, "No ethanolamine or glycerol head group found"

    # Identify fatty acid chain
    # Looking for a chain of at least 20 carbons with multiple double bonds
    fatty_acid_chains = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE and bond.IsInRing() == False:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                # Double bond between carbon atoms
                fatty_acid_chains.append((begin_atom.GetIdx(), end_atom.GetIdx()))

    if len(fatty_acid_chains) < 2:
        return False, "Fatty acid chain does not have enough double bonds"

    # Count total carbons in the longest aliphatic chain
    mol = Chem.AddHs(mol)
    chains = Chem.GetSymmSSSR(mol)
    max_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            paths = Chem.rdmolops.GetShortestPath(mol, atom.GetIdx(), atom.GetIdx())
            if len(paths) > max_chain_length:
                max_chain_length = len(paths)

    if max_chain_length < 20:
        return False, "Fatty acid chain is too short (less than 20 carbons)"

    return True, "Molecule is an endocannabinoid with appropriate fatty acid chain and head group"