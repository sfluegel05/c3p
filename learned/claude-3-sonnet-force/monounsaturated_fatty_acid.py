"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: monounsaturated fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid has:
    - A carboxylic acid group (-C(=O)OH)
    - One double or triple bond in the fatty acid chain
    - Singly bonded carbon atoms in the rest of the chain
    - A continuous hydrocarbon chain of at least 4 carbon atoms

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(carboxyl_pattern) == False:
        return False, "No carboxylic acid group found"

    # Check for one double or triple bond
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    num_triple_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.TRIPLE)
    if num_double_bonds + num_triple_bonds != 1:
        return False, "Not exactly one double or triple bond"

    # Check for hydrocarbon chain
    hydrocarbon_chain = False
    for chain in AllChem.FindAllSubgraphIsomers(mol):
        chain_smiles = Chem.MolToSmiles(chain)
        if "C=C" in chain_smiles or "C#C" in chain_smiles:
            hydrocarbon_chain = True
            break

    if not hydrocarbon_chain:
        return False, "No continuous hydrocarbon chain found"

    # Calculate chain length
    chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            neighbors = atom.GetNeighbors()
            longest_chain = max([len(Chem.FindAllSubgraphIsomers(mol, atomPredicate=AllChem.IsSingleCarbonAtom, startAtom=neighbor)) for neighbor in neighbors], default=0)
            chain_length = max(chain_length, longest_chain)

    if chain_length < 4:
        return False, "Fatty acid chain too short"

    return True, "Contains a carboxylic acid group, one double or triple bond, and a hydrocarbon chain of at least 4 carbon atoms"