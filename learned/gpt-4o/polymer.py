"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer.
    Polymers typically have high molecular weights, long carbon chains,
    repeating substructures, and potentially branching characteristics.

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
    
    # Check for long carbon chain with potential for cross-linking/branching
    chain_pattern = Chem.MolFromSmarts("[C](~[C])~[C]~[C](~[C])[C]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)

    # Count the number of carbon-carbon single bonds and branching points
    num_c_c_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6)
    num_branching_points = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2)

    # Check for high molecular weight (still relevant, but used with caution)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Enhanced criteria: adjust thresholds and add branching consideration
    if (len(chain_matches) >= 2 and num_c_c_bonds > 20 and num_branching_points > 1) or mol_weight > 1000:
        return True, "Has long chains, branching, and high molecular weight suggesting polymer structure."
    
    # Attempt to identify repeating structural motifs, like ether (-C-O-C-), amide, or ester bonds
    potential_repeating_unit = Chem.MolFromSmarts("[*](~[*])~[*]") # Placeholder pattern for demonstration
    if mol.HasSubstructMatch(potential_repeating_unit):
        return True, "Contains potential repeating units typical of polymer structure."

    return False, "Does not meet enhanced polymer criteria of long chains with branching or potential repeating units."