"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI:51638 3-hydroxy fatty acid

Any fatty acid with a hydroxy functional group in the beta- or 3-position. 
beta-Hydroxy fatty acids accumulate during cardiac hypoxia, and can also be used 
as chemical markers of bacterial endotoxins.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondStereo

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for hydroxy group pattern at 3-position
    hydroxy_pattern = Chem.MolFromSmarts("[CX4;H2][CX4;H1][OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    acid_match = mol.GetSubstructMatch(acid_pattern)
    
    if not hydroxy_matches or acid_match == ():
        return False, "Missing hydroxy group at 3-position or carboxylic acid group"
    
    # Get the main aliphatic chain
    chains = rdMolDescriptors.GetLipidChainDescriptors(mol)
    if not chains:
        return False, "No aliphatic chain found"
    main_chain = max(chains, key=lambda c: c.length)
    
    # Check if hydroxy group is attached to main chain
    for match in hydroxy_matches:
        hydroxy_atom = mol.GetAtomWithIdx(match[2])
        hydroxy_bond = hydroxy_atom.GetBondWithIdxIter(match[1]).next()
        if hydroxy_bond.GetStereo() != BondStereo.STEREONONE:
            continue  # Skip if stereochemistry is specified
        
        bond_atoms = [hydroxy_bond.GetBeginAtomIdx(), hydroxy_bond.GetEndAtomIdx()]
        if match[1] in bond_atoms and any(atom.IsInRingSize(3) for atom in mol.GetAtomWithIdx(match[0:2])):
            continue  # Skip if part of a 3-membered ring
        
        if main_chain.startIdx in bond_atoms or main_chain.endIdx in bond_atoms:
            break
    else:
        return False, "Hydroxy group not attached to main aliphatic chain"
    
    return True, "Molecule contains a hydroxy group at the 3-position of the main aliphatic chain"