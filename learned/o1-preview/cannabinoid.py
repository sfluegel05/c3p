"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are a diverse group of compounds that may contain heterocyclic rings with oxygen,
    long hydrocarbon chains, and functional groups such as hydroxyls, amides, esters, or ethers.
    They can be classified into phytocannabinoids, endocannabinoids, and synthetic cannabinoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns
    patterns = []

    # Pattern for heterocyclic ring containing oxygen (e.g., benzopyran ring)
    benzopyran_pattern = Chem.MolFromSmarts("c1cc2ccccc2oc1")
    patterns.append((benzopyran_pattern, "Contains benzopyran ring system"))

    # Pattern for amide linkage to ethanolamine (endocannabinoids like anandamide)
    amide_ethanolamine_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    patterns.append((amide_ethanolamine_pattern, "Contains amide linkage to ethanolamine"))

    # Pattern for glycerol ester of fatty acid (e.g., 2-arachidonoylglycerol)
    glycerol_ester_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)CC")
    patterns.append((glycerol_ester_pattern, "Contains glycerol ester of fatty acid"))

    # Pattern for indole core (synthetic cannabinoids)
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]cc2")
    patterns.append((indole_pattern, "Contains indole core"))

    # Pattern for long hydrocarbon chain (>=10 carbons)
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")  # 10 carbons
    patterns.append((long_chain_pattern, "Contains long hydrocarbon chain"))

    # Pattern for phenolic ring with hydroxyl group
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    patterns.append((phenol_pattern, "Contains phenolic ring with hydroxyl group"))

    # Check each pattern
    for pattern, reason in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, reason

    # As a catch-all, check for any heterocyclic ring containing oxygen
    ring_oxygen = False
    for ring in Chem.GetSymmSSSR(mol):
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        atom_syms = [atom.GetSymbol() for atom in atoms_in_ring]
        if 'O' in atom_syms:
            ring_oxygen = True
            break
    if ring_oxygen:
        # Check for long hydrocarbon chain
        if mol.HasSubstructMatch(long_chain_pattern):
            return True, "Contains heterocyclic ring with oxygen and long hydrocarbon chain"

    return False, "Does not contain structural features characteristic of cannabinoids"