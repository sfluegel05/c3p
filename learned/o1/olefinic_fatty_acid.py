"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is characterized by a long hydrocarbon chain (usually ≥8 carbons)
    containing at least one carbon-carbon double bond (C=C), with a terminal carboxylic acid
    or ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules containing rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not typical for fatty acids"
    
    # Find all acyclic chains in the molecule
    def get_longest_chain(mol):
        longest = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon atom
                paths = Chem.FindAllPathsOfLengthN(mol, mol.GetNumAtoms(), useBonds=True, useHs=False)
                for path in paths:
                    # Check if path is acyclic
                    if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path):
                        if len(path) > len(longest):
                            longest = path
        return longest

    longest_chain = get_longest_chain(mol)
    chain_length = len(longest_chain)
    
    if chain_length < 8:
        return False, f"Longest carbon chain is {chain_length} carbons, requires at least 8"
    
    # Check for at least one carbon-carbon double bond in the longest chain
    has_double_bond = False
    for i in range(len(longest_chain)-1):
        bond = mol.GetBondBetweenAtoms(longest_chain[i], longest_chain[i+1])
        if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
            has_double_bond = True
            break
    if not has_double_bond:
        return False, "No C=C double bond found in the longest carbon chain"
    
    # Check for terminal carboxylic acid or ester group
    # Terminal carboxylic acid group: O=C[O;H1]
    # Ester group: O=C-O-C
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    ester = Chem.MolFromSmarts("C(=O)O[C]")
    if not mol.HasSubstructMatch(carboxylic_acid) and not mol.HasSubstructMatch(ester):
        return False, "No terminal carboxylic acid or ester group found"
    
    return True, "Molecule is an olefinic fatty acid with a long hydrocarbon chain and at least one C=C double bond"