"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer.
    Polymers typically have high molecular weights, long carbon chains,
    and/or repeating substructures.

    Args:
        smiles (str): SMILES string of the chemical entity

    Returns:
        bool: True if the entity is classified as a polymer, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for long carbon chain as a basic criterion
    chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[C]~[C]")  # Representative long carbon chain
    chain_matches = mol.GetSubstructMatches(chain_pattern)

    # Count the number of carbon-carbon single bonds
    num_c_c_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6)
    
    # Establish an arbitrary threshold for considering it as a polymer
    # Here, we give primary consideration to instances of long C chains
    if len(chain_matches) >= 1 and num_c_c_bonds > 10:
        return True, "Has long carbon chains and likely repeating units typical of polymers"
    
    # Check for high molecular weight as polymers tend to have
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    
    if mol_weight > 500:  # Arbitrary threshold for considering high molecular weight 
        return True, "High molecular weight consistent with polymers"

    # Check for potential repeating units
    # A complex task that would involve more robust strategy than simple subgraph searches
    # For simplicity, we'll omit deep searching for repeating units in this basic classification

    return False, "Does not meet basic polymer criteria of long chains or high molecular weight"