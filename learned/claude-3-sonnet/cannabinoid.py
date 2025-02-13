"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:36311 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by the presence of an oxygen atom as part of a heterocyclic ring
    or in various functional groups, long aliphatic chains (typically 15-25 carbon atoms), and one
    or more double bonds in the aliphatic chains.

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

    # Check for presence of oxygen
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8) == 0:
        return False, "No oxygen atom present"

    # Check for long aliphatic chains (15-25 carbon atoms)
    longest_chain = AllChem.GetLongestAliphaticChain(mol)
    if longest_chain is None or (longest_chain.Length() < 15 or longest_chain.Length() > 25):
        return False, "Aliphatic chain length not in the range of 15-25 carbon atoms"

    # Check for presence of double bonds in the aliphatic chain
    if AllChem.CalcNumAliphaticDoubleBonds(mol) == 0:
        return False, "No double bonds present in the aliphatic chain"

    # Check for presence of heterocyclic rings
    heterocyclic_rings = [ring for ring in AllChem.GetSymmSSSR(mol) if ring.IsHeterocyclic()]
    if len(heterocyclic_rings) == 0:
        return True, "Contains oxygen, long aliphatic chain with double bonds (non-cyclic cannabinoid)"

    return True, "Contains oxygen, long aliphatic chain with double bonds, and heterocyclic ring(s)"