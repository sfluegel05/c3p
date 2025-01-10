"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, 'Invalid SMILES string'

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1,H0-]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for epoxide ring (three-membered cyclic ether)
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring found"

    # Check that the epoxide ring is part of an aliphatic chain (not aromatic)
    is_aliphatic_epoxide = False
    for match in epoxide_matches:
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in match]
        if all(not atom.IsInRingOfSize(6) and not atom.GetIsAromatic() for atom in atoms_in_ring):
            is_aliphatic_epoxide = True
            break
    if not is_aliphatic_epoxide:
        return False, "Epoxide ring is not part of an aliphatic chain"

    # Determine the longest aliphatic carbon chain
    from rdkit.Chem.rdmolops import GetMolFrags
    chains = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    max_chain_length = 0
    for chain in chains:
        # Count the number of carbon atoms in the chain
        carbon_count = sum(1 for atom in chain.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() <= 4)
        if carbon_count > max_chain_length:
            max_chain_length = carbon_count

    if max_chain_length < 8:
        return False, f"Longest carbon chain length is {max_chain_length}, which is less than 8 carbons"

    # Ensure that the epoxide ring is connected to the carbon chain
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    epoxide_atoms = set(idx for match in epoxide_matches for idx in match)
    for carboxyl_match in carboxylic_matches:
        carboxyl_carbon_idx = carboxyl_match[0]
        # Use a bond-based search to find if the epoxide is connected
        paths = Chem.rdmolops.GetShortestPath(mol, carboxyl_carbon_idx, list(epoxide_atoms)[0])
        if paths:
            # Check if path consists only of carbons and epoxide oxygens
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() in [6,8] for idx in paths):
                return True, "Molecule contains a carboxylic acid group and an epoxide ring within a hydrocarbon chain"

    return False, "Epoxide ring is not connected to the hydrocarbon chain with carboxylic acid group"