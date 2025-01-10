"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of the 
    carboxy group of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ester linkage pattern: -C(=O)O-
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C,c]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Check for long carbon chains, indicating fatty acid and alcohol
    # A rough estimate: Collect chains attached to both sides of the ester group
    carboxylic_part = Chem.MolFromSmarts("[OX1]=[CX3]([#6])O")
    alcohol_part = Chem.MolFromSmarts("O[#6]")

    carboxy_chain = mol.GetSubstructMatches(carboxylic_part)
    alcohol_chain = mol.GetSubstructMatches(alcohol_part)

    if not carboxy_chain or not alcohol_chain:
        return False, "Insufficient chain length for fatty acid or alcohol"

    # Length of hydrocarbon chains (ignoring functional groups)
    # This logic could be more specific if the molecular representation is complex
    carbon_chain_length = lambda substruct: sum(1 for match in substruct for atom in match if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6)

    # Assume each half needs at least 8 carbons to be considered fatty
    if carbon_chain_length(carboxy_chain) < 8 or carbon_chain_length(alcohol_chain) < 8:
        return False, "Carbon chains are too short to be considered fatty acid/alcohol"

    return True, "Molecule contains a fatty acid ester linkage with sufficient carbon chain length"