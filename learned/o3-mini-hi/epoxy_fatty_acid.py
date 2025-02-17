"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    The molecule must contain a carboxylic acid group (fatty acid) and an epoxide ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an epoxy fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a carboxylic acid group.
    # This pattern will match a carbon with a double-bond to oxygen and a single bond to an -OH.
    ca_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(ca_pattern):
        return False, "No carboxylic acid group found to signify a fatty acid"
    
    # Define the SMARTS for an epoxide ring:
    # Look for a three-membered ring consisting of two carbons and one oxygen.
    epoxide_pattern = Chem.MolFromSmarts("[C;R3][O;R3][C;R3]")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"
    
    # Count the total number of carbon atoms in the molecule.
    # Fatty acids typically have long carbon chains.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, f"Too few carbon atoms ({c_count}) to be a typical fatty acid"

    # (Optional) Check approximate molecular weight. Many fatty acids are in the range above 200 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.2f} Da) too low for a fatty acid"
    
    # If all checks pass, we classify the molecule as an epoxy fatty acid.
    return True, "Contains a carboxylic acid group with a sufficiently long aliphatic chain and an epoxide ring indicative of an epoxy fatty acid"

# Example usage (uncomment to test):
# smiles_example = "CCCCCCCCCCCC(=O)O[C@H]1[C@@H](O)[C@H](C)O1"  # (Example; adjust to valid epoxy fatty acid SMILES)
# result, reason = is_epoxy_fatty_acid(smiles_example)
# print(result, reason)